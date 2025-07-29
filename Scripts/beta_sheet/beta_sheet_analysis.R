cat("Running: beta_sheet_analysis.R\n")

############################################################################################

# Installing and Loading Packages
required_packages <- c("dplyr", "tidyr","readr", "ggplot2", "bio3d", "zoo", "scales")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]

# Install missing packages
if(length(not_installed)) install.packages(not_installed)

# Load all required packages
invisible(lapply(required_packages, require, character.only = TRUE))

# Cleaning Environment
rm(list = ls())

############################################################################################

# Creating save dir

if(!dir.exists(file.path("Output", "b_sheets"))) dir.create(file.path("Output", "b_sheets"))
if(!dir.exists(file.path("Output", "b_sheets", "ramachandran_plots"))) dir.create(file.path("Output", "b_sheets", "ramachandran_plots"))

############################################################################################

# Specifying range for B-sheet conformation
# References: 
# B-Turns and their distortions: a proposed new nomenclature - Figure 2
# A fresh look at the Ramachandran plot and the occurrence of standard structures in proteins (used to find above paper)
psi_min = 60
psi_max = 180
psi_min_extended = -180
psi_max_extended = -150
phi_min = -180
phi_max = -90

# Initialising dataframe
all_torsion_data <- data.frame()
proportion_data <- data.frame()

# Locate all pdbs in folder that need to be looped through
pdb_path <- file.path("Output", "PDBs", "published_structure")
# Get a list of all .pdb files in the folder
pdb_files <- list.files(path = pdb_path, pattern = "\\.pdb$", full.names = TRUE)

### Removing PDBs with a poor mean resolution - identified in the validation scripts ###

# Getting the names of the PDB files in the asymmetric units folder
pdb_filenames <- basename(pdb_files)
pdb_filenames <- sub("\\.pdb$", "", pdb_filenames)
pdb_names <- sub("_.*", "", pdb_filenames)

# Getting the list of high resolution PDBs
high_resolution_PDBs <- read_delim(file.path("Output", "Validation", "high_resolution_pdb_ids.txt"), 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
# Creating PDB column
high_resolution_PDBs$pdb <- sub("_.*", "", high_resolution_PDBs$pdb_id)


# Identifying the index of PDBs in pdb_names that are not found in high_resolution_PDBs
low_res <- setdiff(pdb_names, high_resolution_PDBs$pdb)
# Finding the indices of elements in pdb_names that match the values in low_res
indices_to_remove <- which(pdb_names %in% low_res)
# Removing the corresponding elements from pdb_files
pdb_files_filtered <- pdb_files[-indices_to_remove]

for (i in 1:length(pdb_files_filtered)) {
  
  cat("Analysing -", pdb_names[i], "\n")  
  
  # Reading in PDB file
  pdb <- read.pdb(pdb_files_filtered[i])
  
  # Check if the PDB file contains any atoms
  if (!is.null(pdb$atom) && nrow(pdb$atom) > 0) {
   
    # Get all of the torsion angles
    torsions <- torsion.pdb(pdb)
    
    # Convert the torsions list into a dataframe
    torsion_df <- as.data.frame(torsions$tbl)
    
    # Residue Number, Chain and residue type are stored in the row names so making into new columns
    torsion_df$info <- rownames(torsion_df)
    torsion_df <- separate(torsion_df, info, into = c("pos", "chain", "residue"), sep = "\\.")
    torsion_df$pos <- as.numeric(torsion_df$pos)
    
    # Adding in PDB name
    torsion_df$pdb <- pdb_names[i]
    
    # Removing NAs and UNK residues
    torsion_df <- torsion_df %>% filter(residue != "UNK")
    torsion_df <- torsion_df %>% filter(!(is.na(psi) & is.na(phi)))
    
    # Adding to all_torsion_data
    all_torsion_data <- rbind(all_torsion_data, torsion_df)
    
    # Ramachandran Plot
    ggplot(torsion_df, aes(x = phi, y = psi)) +
      annotate("rect", xmin = phi_min, xmax = phi_max, ymin = psi_min, ymax = psi_max, fill = "grey", alpha = 0.75) +
      annotate("rect", xmin = phi_min, xmax = phi_max, ymin = psi_min_extended, ymax = psi_max_extended, fill = "grey", alpha = 0.75) +
      annotate("segment", x = -180, xend = 180, y = 0, yend = 0, linewidth = 0.75) +
      annotate("segment", x = 0, xend = 0, y = -180, yend = 180, linewidth = 0.75) +
      geom_point(color = "blue", alpha = 0.7, size = 4) +
      scale_x_continuous(breaks = seq(-180, 180, 30), limits = c(-180, 180), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(-180, 180, 30), limits = c(-180, 180), expand = c(0, 0)) +
      theme_bw() +
      theme(plot.title = element_text(size = 24, colour = "black", face = "bold", hjust = 0.5),
            axis.text = element_text(size = 16, colour = "black"),
            axis.title = element_text(size = 20, colour = "black", face = "bold")) +
      labs(title = paste0("Ramachandran Plot - ", pdb_names[i]),
           x = "Phi (φ) Angles", y = "Psi (ψ) Angles")

    # Saving plots
    ggsave(plot = last_plot(), file = paste0("Output/b_sheets/ramachandran_plots/", pdb_names[i], ".png"),
           height = 8, width = 8)

    # Filter for residues by torsion angles
    b_strands <- torsion_df[torsion_df$psi > psi_min & 
                              torsion_df$psi < psi_max & 
                              torsion_df$phi > phi_min & 
                              torsion_df$phi < phi_max |
                              torsion_df$psi > psi_min_extended & 
                              torsion_df$psi < psi_max_extended & 
                              torsion_df$phi > phi_min & 
                              torsion_df$phi < phi_max, ]
    
    ### Calculating the proportion of molecules a residue is found in a b-strand ###
    
    # Counting unique chains
    num_chains <- as.numeric(length(unique(torsion_df$chain)))
    
    # Counting occurences of each residue
    residue_counts <- b_strands %>% count(pos)
    
    # Removing NA's
    residue_counts <- na.omit(residue_counts)
    
    # Calculating the proportion of monomers a residue is found to be in a beta sheet
    residue_counts$num_chains <- num_chains
    residue_counts$proportion <- residue_counts$n / residue_counts$num_chains
    
    # Adding in reference column
    residue_counts$pdb <- pdb_names[i]
    
    # Merging current data to overall data frame
    proportion_data <- rbind(proportion_data, residue_counts)
    
    
  } else {
    cat("Skipping ", pdb_names[i], ": PDB file is empty or contains no atomic coordinates.")
  }

}

# Saving data
write.csv(all_torsion_data, file = file.path("Output", "b_sheets", "all_torsion_data.csv"), row.names = FALSE)
write.csv(proportion_data, file = file.path("Output", "b_sheets", "proportion_in_b_strand.csv"), row.names = FALSE)

### Identify 4+ b-strands ###
find_consecutive <- function(vec, min_run = 4) {
  rle_vals <- rle(diff(vec) == 1)  # Find sequences where the difference is 1
  run_lengths <- rep(rle_vals$lengths + 1, rle_vals$lengths)  # Expand run lengths
  
  # Mark values as part of a run if they belong to a sequence of `min_run` or more
  in_run <- run_lengths >= min_run
  c(FALSE, in_run) | c(in_run, FALSE)  # Extend check to cover full run
}

# Sort the data (just in case) and apply the function
proportion_data_four_plus <- proportion_data[find_consecutive(proportion_data$pos), ]

# Saving data
write.csv(proportion_data_four_plus, file = file.path("Output", "b_sheets", "proportion_in_b_strand_four_plus.csv"), row.names = FALSE)

#######################################################
########################################################

# Ramachandran Plot for all pdbs to see overall spread
ggplot(all_torsion_data, aes(x = phi, y = psi)) +
  annotate("rect", xmin = phi_min, xmax = phi_max, ymin = psi_min, ymax = psi_max, fill = "black", alpha = 0.5) +
  annotate("rect", xmin = phi_min, xmax = phi_max, ymin = psi_min_extended, ymax = psi_max_extended, fill = "black", alpha = 0.5) +
  annotate("segment", x = -180, xend = 180, y = 0, yend = 0, linewidth = 0.9) +
  annotate("segment", x = 0, xend = 0, y = -180, yend = 180, linewidth = 0.9) +
  geom_point(color = "blue", alpha = 0.2, size = 2) +
  scale_x_continuous(breaks = seq(-180, 180, 30), limits = c(-180, 180), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-180, 180, 30), limits = c(-180, 180), expand = c(0, 0)) +
  theme_bw() +
  theme(plot.title = element_text(size = 24, colour = "black", face = "bold", hjust = 0.5),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 20, colour = "black", face = "bold")) +
  labs(title = paste0("Ramachandran Plot"), 
       x = "Phi (φ) Angles", y = "Psi (ψ) Angles") 

ggsave(plot = last_plot(), filename = file.path("Output", "b_sheets", "all_pdbs_ramachandran.png"),
       height = 8, width = 8)

# Removing duplicate residues grouped by PDB to account for different numbers of published chains
unique_all_torsion_data <- all_torsion_data %>%
                             group_by(pdb) %>%
                              filter(!duplicated(pos)) %>%
                                ungroup()

# Counting the number of occurrences of each residue
total_counts <- unique_all_torsion_data %>% count(pos)
names(total_counts)[names(total_counts) == "n"] <- "total"

#################################################################
### Calculating the total b-strand proportion across all pdbs ###
#################################################################

total_PDBs <- length(pdb_files_filtered)

# 1+
total_b_strand <- proportion_data %>% 
                  group_by(pos) %>%
                  summarise(num_b_strand = sum(proportion))

total_b_strand$percentage <- (total_b_strand$num_b_strand / total_PDBs) * 100

# 4+
total_b_strand_four_plus <- proportion_data_four_plus %>% 
                            group_by(pos) %>%
                            summarise(num_b_strand_four_plus = sum(proportion))

total_b_strand_four_plus$percentage_four_plus <- (total_b_strand_four_plus$num_b_strand_four_plus / total_PDBs) * 100


# Merging total, 1+ and 4+ counts and calculating the proportion for each
plot_data <- merge(total_counts, total_b_strand, by = "pos")
plot_data$percentage <- (plot_data$num_b_strand / plot_data$total) * 100
plot_data <- merge(plot_data, total_b_strand_four_plus, by = "pos")
plot_data$percentage_four_plus <- (plot_data$num_b_strand_four_plus / plot_data$total) * 100

### Errico et al 2025 used a sliding window of 7 so repeating to compare ###
### https://doi.org/10.1042/bcj20240602 ### 

# Sliding window of percentage
window_size = 3

plot_data$sliding_window <- rollapply(plot_data$percentage, width = window_size, FUN = mean, align = "center", fill = NA)
plot_data$sliding_window_four_plus <- rollapply(plot_data$percentage_four_plus, width = window_size, FUN = mean, align = "center", fill = NA)

# Plotting counts
ggplot(data = plot_data) +
  geom_line(aes(x = pos, y = percentage, linetype = "1+"), linewidth = 0.75) +
  geom_line(aes(x = pos, y = percentage_four_plus, linetype = "4+"), linewidth = 0.75, colour = "red") +
  ylab(bquote(bold("Fraction in " * beta * "-conformation (%)"))) +
  xlab("Residue Number") +
  scale_x_continuous(breaks = seq(0,999,20)) +
  scale_y_continuous(breaks = seq(0,100,20)) +
  scale_linetype_manual(values = c("1+" = "solid", 
                                   "4+" = "solid"),
                        label = c(bquote(bold(L[beta] * " \u2265 " * 1)),
                                  bquote(bold(L[beta] * " \u2265 " * 4)))) +
  guides(linetype = guide_legend(override.aes = list(size = 6))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.position = "right")

ggsave(plot = last_plot(), filename = file.path("Output", "b_sheets", "strands.png"),
       height = 4, width = 10)

# Plotting counts
ggplot(data = plot_data) +
  geom_line(aes(x = pos, y = sliding_window, linetype = "1+"), linewidth = 0.75) +
  geom_line(aes(x = pos, y = sliding_window_four_plus, linetype = "4+"), linewidth = 0.75, colour = "red") +
  ylab(bquote(bold("Fraction in " * beta * "-conformation (%)"))) +
  xlab("Residue Number") +
  scale_x_continuous(breaks = seq(0,999,20)) +
  scale_y_continuous(breaks = seq(0,100,20)) +
  scale_linetype_manual(values = c("1+" = "solid", 
                                   "4+" = "solid"),
                        label = c(bquote(bold(L[beta] * " \u2265 " * 1)),
                                  bquote(bold(L[beta] * " \u2265 " * 4)))) +
  guides(linetype = guide_legend(override.aes = list(size = 6))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.position = "right")

ggsave(plot = last_plot(), filename = file.path("Output", "b_sheets", paste0("strands_sliding_window_", window_size, ".png")),
       height = 4, width = 10)

# Importing stablising regions
stable_regions <- read_csv(file.path("Output", "stable_regions", "stable_regions.csv"))

# Splitting region into min and max residue number for plotting reasons
stable_regions <- separate(stable_regions, residues, into = c("min", "max"), sep = "-")

stable_regions$min <- as.numeric(stable_regions$min)
stable_regions$max <- as.numeric(stable_regions$max)
stable_regions$region <- as.character(stable_regions$region)

stable_regions$region <- factor(stable_regions$region, levels = sort(as.numeric(unique(stable_regions$region))))

# Assigning colours
colours <- c("red", "blue", "gold2", "purple", "darkorange2", "cyan", "green", "pink", "chocolate4", 
             "darkgreen", "magenta", "deepskyblue", "darkviolet", "coral", "darkslategray")

# Extend the color vector to match `region`
region_colours <- rep(colours, length.out = nrow(stable_regions))

### Comparing B-sheet propensity and foldx deltaG ###

foldx_stability <- read_csv(file.path("Output", "thermodynamics", "foldx_stability.csv"))

# Finding the average deltaG for each residue position across all PDBs
foldx_df <- foldx_stability %>% 
              group_by(Pos) %>%
                summarise(mean_energy = mean(total, na.rm = TRUE))  

# Inverting DeltaG so increasing stability is in the same direction as increasing % B-sheet
foldx_df$inverted_mean_energy <- foldx_df$mean_energy * -1

# Calculating sliding window
foldx_df$inverted_sliding_window <- rollapply(foldx_df$inverted_mean_energy, width = window_size, FUN = mean, align = "center", fill = NA)

# Scaling B-sheet and DeltaG so they can be compared
foldx_df$scaled_inverted_sliding_window <- rescale(foldx_df$inverted_sliding_window)
plot_data$scaled_sliding_window <- rescale(plot_data$sliding_window)
plot_data$scaled_sliding_window_four_plus <- rescale(plot_data$sliding_window_four_plus)

# Plotting
ggplot() +
  geom_rect(data = stable_regions, aes(xmin = min, xmax = max, 
                                       ymin = -Inf, ymax = Inf, 
                                       fill = region), alpha = 0.75) +
  geom_line(data = plot_data, aes(x = pos, y = scaled_sliding_window, linetype = "B-Strands"), 
            size = 0.75) +
  geom_line(data = foldx_df, aes(x = Pos, y = scaled_inverted_sliding_window, linetype = "DeltaG"), 
            size = 0.75) +
  scale_x_continuous(breaks = seq(0,9999,20)) +
  scale_y_continuous(
    name = bquote(bold("Scaled Fraction in " * beta * "-conformation")),
    sec.axis = sec_axis(~ ., name = bquote(bold("Inverted Scaled Mean " * Delta * "G"*degree * " (kcal.mol"^-1 * ")")))) +
  scale_fill_manual(values = region_colours, guide = "none") +
  scale_linetype_manual(values = c("B-Strands" = "dashed", 
                                   "DeltaG" = "solid"),
                        name = "Data",
                        label = c(bquote(bold(L[beta])),
                                  bquote(bold(Delta * "G"*degree)))) + 
  guides(linetype = guide_legend(override.aes = list(size = 8))) +
  ggtitle(paste0("Sliding Window of ", window_size)) +
  xlab("Residue Number") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        axis.text = element_text(size = 10, face = "bold", colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.position = "right") 

ggsave(plot = last_plot(), filename = file.path("Output", "b_sheets", paste0("strands_vs_deltaG_", window_size, ".png")),
       height = 4, width = 10)

# Plotting
ggplot() +
  geom_rect(data = stable_regions, aes(xmin = min, xmax = max, 
                                       ymin = -Inf, ymax = Inf, 
                                       fill = region), alpha = 0.75) +
  geom_line(data = plot_data, aes(x = pos, y = scaled_sliding_window_four_plus, linetype = "B-Strands"), 
            size = 0.75) +
  geom_line(data = foldx_df, aes(x = Pos, y = scaled_inverted_sliding_window, linetype = "DeltaG"), 
            size = 0.75) +
  scale_x_continuous(breaks = seq(0,9999,20)) +
  scale_y_continuous(
    name = bquote(bold("Scaled Fraction in " * beta * "-conformation")),
    sec.axis = sec_axis(~ ., name = bquote(bold("Inverted Scaled Mean " * Delta * "G"*degree * " (kcal.mol"^-1 * ")")))) +
  scale_fill_manual(values = region_colours, guide = "none") +
  scale_linetype_manual(values = c("B-Strands" = "dashed", 
                                   "DeltaG" = "solid"),
                        name = "Data",
                        label = c(bquote(bold(L[beta] * " \u2265 " * 4)),
                                bquote(bold(Delta * "G"*degree)))) + 
  guides(linetype = guide_legend(override.aes = list(size = 8))) +
  ggtitle(paste0("Sliding Window of ", window_size)) +
  xlab("Residue Number") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        axis.text = element_text(size = 10, face = "bold", colour = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.position = "right") 

ggsave(plot = last_plot(), filename = file.path("Output", "b_sheets", paste0("strands4plus_vs_deltaG_", window_size, ".png")),
       height = 4, width = 10)

# Saving combined plot data
names(foldx_df)[names(foldx_df) == "Pos"] <- "pos"

thermo_and_b_sheet <- merge(plot_data, foldx_df, by = "pos", all = TRUE)

# Selecting columns
thermo_and_b_sheet <- thermo_and_b_sheet %>% select(c("pos", "scaled_sliding_window", "scaled_sliding_window_four_plus", "scaled_inverted_sliding_window"))

names(thermo_and_b_sheet) <- c("Residue", "proportion_in_b_conformation", "proportion_in_b_conformation_four_plus", "scaled_inverted_deltaG")

# Saving thermo_and_b_sheet
write.csv(thermo_and_b_sheet, file.path("Output", "b_sheets", paste0("thermo_vs_b-sheets_window_", window_size, ".csv")), row.names = FALSE)

