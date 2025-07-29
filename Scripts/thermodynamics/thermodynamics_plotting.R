cat("\nRunning: thermodynamics_plotting.R\n")

# Installing and Loading Packages
required_packages <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "zoo", "ggrepel", "dendextend", "FSA")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]   

# Install missing packages
if(length(not_installed)) install.packages(not_installed)    

# Load all required packages
invisible(lapply(required_packages, require, character.only = TRUE))

# Cleaning Environment
rm(list = ls())

################################################################################

### Reading in configuration from run_analysis.sh ###

# Function to read configuration from stdin
read_config <- function() {
  config <- list()
  lines <- readLines(file("stdin"))
  for (line in lines) {
    parts <- strsplit(line, "=")[[1]]
    config[[parts[1]]] <- parts[2]
  }
  return(config)
}

# Read the configuration from stdin
config <- read_config()

remove_poorly_resolved <- as.numeric(config$remove_poorly_resolved)

################################################################################

# Creating Directories
if(!dir.exists(file.path("Output", "thermodynamics", "deltaG_analysis"))) dir.create(file.path("Output", "thermodynamics", "deltaG_analysis"))
if(!dir.exists(file.path("Output", "thermodynamics", "single_PDB_plots"))) dir.create(file.path("Output", "thermodynamics", "single_PDB_plots"))
if(!dir.exists(file.path("Output", "thermodynamics", "residue_analysis"))) dir.create(file.path("Output", "thermodynamics", "residue_analysis"))
if(!dir.exists(file.path("Output", "thermodynamics", "pearson_correlation"))) dir.create(file.path("Output", "thermodynamics", "pearson_correlation"))
if(!dir.exists(file.path("Output", "thermodynamics", "deltaG_comparisons"))) dir.create(file.path("Output", "thermodynamics", "deltaG_comparisons"))

################################################################################

### Importing Data ###

df <- read_csv(file.path("Output", "thermodynamics", "foldx_stability.csv"))

# Removing non-cryoEM structures
metadata <- read_delim(file.path("Output", "selected_pdbs_metadata.txt"), 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

# Selecting results that did not use cryoEM
tmp <- metadata %>% filter(Method != "cryoEM") 

# Filtering df to remove non-cryoEM structures
df <- df %>% filter(!(PDB %in% tmp$`PDB ID`))

################################################################################

### Removing regions outside the fibril core ###
metadata <- subset(metadata, select = c("PDB ID",  "Residues Ordered"))
names(metadata) <- c("PDB", "ordered_residues")

metadata$core_start <- NA
metadata$core_end <- NA


setdiff(metadata$PDB, df$PDB)


for (i in 1:nrow(metadata)) {
  
  range <- metadata$ordered_residues[i]
  range <- str_split_1(range, paste("-", ",", sep = "|"))
  
  metadata$core_start[i] <- as.numeric(range[1])
  metadata$core_end[i] <- as.numeric(range[length(range)])
  
}

# Combining df and metadata
df <- merge(df, metadata, by = "PDB")

# Changing Start and End NAs for non-amyloid atlas pdbs
df <- df %>%
  mutate(core_start = replace_na(core_start, 0),
         core_end = replace_na(core_end, Inf))

# Filtering residues that are outside the core
df <- df %>% filter(Pos >= core_start, 
                    Pos <= core_end)

# Removing columns 
df <- df %>% select(-c("ordered_residues", "core_start", "core_end"))

################################################################################

### Adding in which fibril each chain belongs to ###

# Importing data specifying which fibril each chain is on (generated in extend_fibril_layers.py)
chain_fibril <- read_csv(file.path("Output", "PDBs", "fibrils_extended", "chain_fibril.csv"))

# Renaming columns to match
names(df)[names(df) == "Mol"] <- "chain"

# Creating new fbril column in df based on matches in chain_fibril
df <- left_join(df, chain_fibril, by = c("PDB", "chain"))

################################################################################

### Removing Head and Tail Chains ###

# Importing data showing which chains are head and tails for each fibril
exterior_chains <- read_csv(file.path("Output", "PDBs", "fibrils_extended", "exterior_chains.csv"))

filtered_df <- df[!(paste(df$PDB, df$chain) %in% paste(exterior_chains$PDB, exterior_chains$chain)), ]

################################################################################

### Adding in pdb_id to handle fibrils with 2+ distinct polymorphs ###

# Importing data
COM_and_fibril <- read.csv(file.path("Output", "PDBs", "COM_and_fibril.csv"))

# Selecting columns of interest
pdb_id_data <- COM_and_fibril %>% select(c("PDB", "pdb_id", "fibril"))

# Removing duplicate rows
pdb_id_data <- pdb_id_data %>% distinct()

# Adding pdb_id to filtered_df
filtered_df <- merge(pdb_id_data, filtered_df, by = c("PDB", "fibril"))

################################################################################

if (remove_poorly_resolved == 1) {
  
  ##############################################
  ### Removing residues with poor resolution ###
  ##############################################
  
  # Importing data
  high_resolution_residues <- read.csv(file.path("Output", "Validation", "high_resolution_residues.csv"))
  
  # Renaming columns to match df
  names(high_resolution_residues)[names(high_resolution_residues) == "resno"] <- "Pos"
  
  # Merging with main data frame
  filtered_df <- merge(filtered_df, high_resolution_residues, by = c("pdb_id", "fibril", "Pos"))

}

if (remove_poorly_resolved == 0) {
  
  ##############################################################
  ### Removing only low residue PDBs not individual residues ###
  ##############################################################
  
  # Getting a list of PDBs with high resolution
  high_resolution_PDBs <- unique(high_resolution_residues$pdb_id)
  
  # Keeping only PDBs found in high_resolution_PDBs
  filtered_df <- filtered_df %>% filter(pdb_id %in% high_resolution_PDBs)  
}


################################################################################

# Saving filtered_df
write.csv(filtered_df, file = file.path("Output", "thermodynamics", "foldx_stability_filtered.csv"), row.names = FALSE)

# Cleaning Environment
rm(list = setdiff(ls(), c("filtered_df", "df")))
invisible(gc())


###################################################################################################
###################################################################################################
###################################################################################################

### Total Energy Contribution ###

###################################################################################################
###################################################################################################
###################################################################################################


# Calculating the mean delta G per residue for each PDB 
mean_per_res <- filtered_df %>%
                  group_by(pdb_id, Pos) %>%
                    summarize(mean_energy = mean(total, na.rm = TRUE)) %>%
                      ungroup()

# Saving mean_per_res
write.csv(mean_per_res, file = file.path("Output", "thermodynamics", "mean_deltaG_per_residue.csv"), row.names = FALSE)


# Calculating the threshold for residues to be considered stabilising
mean_stability <- mean(mean_per_res$mean_energy, na.rm = TRUE)
SD_stability <- sd(mean_per_res$mean_energy, na.rm = TRUE)
stabilising_threshold <- mean_stability - SD_stability
destabilising_threshold <- mean_stability + SD_stability


# Plotting

### Density Plot ###
ggplot() +
  geom_density(data = mean_per_res, aes(x = mean_energy, group = pdb_id), 
               colour=alpha("black", 0.25), linewidth = 1) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 1, linetype = "dashed", alpha = 0.9) +
  ylab("Density") +
  xlab(bquote(bold(Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
  scale_x_continuous(breaks = seq(-100,100, by = 1)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text.y = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 1.5), 
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.line.x = element_line(color="black", linewidth = 1.5),
        axis.line.y = element_line(color="black", linewidth = 1.5),
        legend.position = "none",
        legend.title = element_text(size = 20, colour = "black", face = "bold"),
        legend.text = element_text(size = 16, colour = "black", face = "bold")) +
  labs(colour = "PDB")

ggsave(plot = last_plot(), file = file.path("Output", "thermodynamics", "deltaG_analysis", "stability_density_plot.png"),
       height = 6, width = 8)

### Line Plot ###
ggplot(data = mean_per_res, aes(x=Pos)) +
  geom_line(aes(y = mean_energy, group = interaction(pdb_id, cumsum(c(0, diff(Pos) != 1)))), 
            linewidth = 0.75, alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 0.75, linetype = "solid", alpha = 0.9) +
  ylab(bquote(bold(Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
  xlab("Residue Number") +
  scale_x_continuous(breaks = seq(0, 99999, by = 1), 
                     labels = ifelse(seq(0, 99999, by = 1) %% 10 == 0, seq(0, 99999, by = 1), ""),
                     expand = expansion(mult = c(0.01,0.01))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text = element_text(size = 14, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 1.5), 
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.position = "none")

ggsave(plot = last_plot(), file = file.path("Output", "thermodynamics", "deltaG_analysis", "stability_line_plot.png"),
       height = 6, width = 12)


############################################################################################

### Plotting each PDB individually so I can make a movie and visualize the shape of each ###

############################################################################################

pdb_names <- unique(filtered_df$pdb_id)

for (pdb in pdb_names) {

  x <- mean_per_res %>% filter(pdb_id %in% pdb)

  ggplot(data = x, aes(x=Pos)) +
    geom_line(aes(y=mean_energy, group = cumsum(c(0, diff(Pos) != 1))), linewidth = 2) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 1, linetype = "dashed", alpha = 0.9) +
    geom_hline(yintercept = -0.5, colour = "red", linewidth = 1, linetype = "dashed", alpha = 0.9) +
    ggtitle(paste(pdb, "Sliding Window of 5 Mean", sep = " - ")) +
    ylab(bquote(bold(Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
    xlab("Residue Number") +
    #scale_x_continuous(limits = c(0, 100), breaks = seq(0,99999, by = 10)) +
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          panel.grid.major = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
          plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
          axis.text.y = element_text(size = 14, colour = "black", face = "bold"),
          axis.text.x = element_text(size = 14, colour = "black", face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold", vjust = 1.5),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.line.x = element_line(color="black", linewidth = 1.5),
          axis.line.y = element_line(color="black", linewidth = 1.5))
  

  ggsave(plot = last_plot(), file = file.path("Output", "thermodynamics", "single_PDB_plots", paste0(pdb, "_stability.png")),
         height = 6, width = 12)

}

############################################################################################

### Calculating the Pearson Correlation Coefficient Between PDB Stability ### 

############################################################################################ 

wide_df <- mean_per_res %>%
  select(c("pdb_id", "Pos", "mean_energy")) %>%
  pivot_wider(names_from = pdb_id, values_from = mean_energy) 

# Select columns for correlation calculation (excluding the first column)
columns_to_correlate <- wide_df[, -1]

# Calculate correlation matrix, ignoring NA values
correlation_matrix <- cor(columns_to_correlate, use = "pairwise.complete.obs")

# Convert correlation matrix to data frame
correlation_df <- as.data.frame(as.table(correlation_matrix))
names(correlation_df) <- c("Var1", "Var2", "Correlation")

# Plotting Correlation
ggplot(correlation_df, aes(Var1, Var2, fill = Correlation)) +
  geom_tile(linewidth = 0.2, colour = "black") +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       values = scales::rescale(c(-1, 0, 1)),
                       breaks = c(-1,-0.5,0,0.5,1),
                       labels = c(-1,-0.5,0,0.5,1),
                       limits = c(-1,1)) +
  labs(title = "Pearson Correlation Between\nPDBs Per Residue FoldX Stability", x = "", y = "") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", 
                                   angle = 90, vjust = 0.25, hjust = 1),
        axis.line = element_blank(),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        legend.text = element_text(size = 10, colour = "black"))

ggsave(plot = last_plot(), file = file.path("Output", "thermodynamics", "pearson_correlation", "Pearson_heatmap.png"),
       height = 11, width = 12, bg = "white")

write.csv(correlation_df, file = file.path("Output", "thermodynamics", "pearson_correlation", "Pearson_Scores.csv"), row.names = FALSE)


### Correlation Box Plot ###

# Adding a dummy column to use for x-axis
correlation_df$dummy <- ""

# Creating a comparison column
correlation_df <- correlation_df %>%
                    mutate(comparison = paste0(Var1, "-", Var2))

# Removing Self Comparisons
correlation_df <- correlation_df %>% filter(Var1 != Var2)

label <- bquote(bold(atop("Pearson correlation between",  
                          "PDBs " * Delta * G^degree * 
                            " per residue (kcal·mol"^-1 * ")")))

# plot
ggplot(data = correlation_df, aes(x = dummy, y = Correlation)) +
  geom_violin(fill = "grey70", linewidth = 1, width = 0.5) +
  geom_boxplot(size = 1, colour = "black", width = 0.1, fill = NA, outliers = FALSE) +
  scale_y_continuous(name = label, limits = c(-1,1)) +
  xlab("") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.3, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.3, linetype = "dashed"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text = element_text(size = 14, colour = "black", face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        axis.line = element_line(color = "black", linewidth = 1.5),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        legend.text = element_text(size = 10, colour = "black"))

ggsave(plot = last_plot(), file = file.path("Output", "thermodynamics", "pearson_correlation", "Pearson_boxplot.png"),
       height = 6, width = 5)

### Saving Median and Mean Correlation ###
sink(file.path("Output", "thermodynamics", "pearson_correlation", "average_correlation.txt"))
cat("Average Pearson Correlation Score\n\n")
cat("Mean: ", mean(correlation_df$Correlation, na.rm = TRUE), "\n")
cat("Median: ", median(correlation_df$Correlation, na.rm = TRUE), "\n")
cat("SD: ", sd(correlation_df$Correlation, na.rm = TRUE), "\n")
sink()

##################################################################################################
##################################################################################################
##################################################################################################

### Plotting individual energies ### 

##################################################################################################
##################################################################################################
##################################################################################################

rm(list = setdiff(ls(), c("filtered_df", "stabilising_threshold", "destabilising_threshold")))


# Average per PDB and Residue

selected_columns <- c("total", "sideHbond", "energy_VdW", "electro", "energy_SolvP", "energy_SolvH", "energy_vdwclash")

selected_columns <- colnames(filtered_df[11:32])

# Calculating the mean delta G per residue for each PDB 
pdb_means <- filtered_df %>%
             group_by(pdb_id, Pos) %>%
             summarise(across(all_of(selected_columns), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"))


# Reshape the dataframe from wide to long format
long_df <- pdb_means %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "variable", values_to = "value") 


# Removing NA values
long_df <- na.omit(long_df)

# Remove "mean_" pattern
long_df$variable <- gsub("mean_", "", long_df$variable)

# Setting order of long_df (important for cumsum() to calculate line breaks)
long_df <- long_df[order(long_df$variable, long_df$pdb_id, long_df$Pos), ]

# Setting plot order
long_df$variable <- factor(long_df$variable, levels = selected_columns)

# Create the plot
ggplot(long_df, aes(x = Pos, y = value)) +
  geom_line(aes(group = interaction(pdb_id, cumsum(c(0, diff(Pos) != 1)))), 
            linewidth = 1, alpha = 0.25) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  facet_wrap(~variable) +
  #labs(x = "Position", y = "\u0394G per Residue (kcal.mol-1)", 
  #     title = "Mean Scores - Separate PDBs") +
  ggtitle("Mean Scores - Separate PDBs") +
  ylab(bquote(bold(Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
  xlab("Position") +
    theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.line = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

ggsave(plot = last_plot(), file = file.path("Output", "thermodynamics", "deltaG_analysis", "all_energies_separate_pdb.png"),
       height = 12, width = 15, bg = "white")


###################################
### Characterizing Residue Type ### 
###################################

filtered_df$Code[filtered_df$Code == "H1S"] <- "HIS"
filtered_df$Code[filtered_df$Code == "H2S"] <- "HIS"

# Counting how many PDBs each residue position is found in
all_pos_counts <- filtered_df %>% count(Pos)

# Counting how many PDBs each residue position is found to be stabilising
stabilising_df <- filtered_df %>% filter(total <= stabilising_threshold)
destabilising_df <- filtered_df %>% filter(total >= destabilising_threshold)

# Creating a amino acid property data frame
amino_acids <- data.frame(
  Code = c("ALA", "ARG", "ASN", "ASP", 
           "CYS", "GLN", "GLU", "GLY", 
           "HIS", "ILE", "LEU", "LYS", 
           "MET", "PHE", "PRO", "SER", 
           "THR", "TRP", "TYR", "VAL"),
  property = c("Hydrophobic", "Positive", "Polar", "Negative", 
               "Polar", "Polar", "Negative", "Polar", 
               "Positive", "Hydrophobic", "Hydrophobic", "Positive", 
               "Hydrophobic", "Aromatic", "Hydrophobic", "Polar", 
               "Polar", "Aromatic", "Aromatic", "Hydrophobic")
)


# Counting the total of each residue type for all PDBs

all_residue_counts <- filtered_df %>% count(Code)

# Counting the residue type from stabilising residues for all PDBs
stable_residue_counts <- stabilising_df %>% count(Code)

# Counting the residue type from destabilising residues for all PDBs
destable_residue_counts <- destabilising_df %>% count(Code)

# Normalising the stablising counts by the total counts
names(all_residue_counts)[names(all_residue_counts) == "n"] <- "total_count"
names(stable_residue_counts)[names(stable_residue_counts) == "n"] <- "stabilising_count"
names(destable_residue_counts)[names(destable_residue_counts) == "n"] <- "destabilising_count"

# Merging Data
residue_count_df <- merge(stable_residue_counts, destable_residue_counts, by = "Code", all = TRUE)
residue_count_df <- merge(residue_count_df, all_residue_counts, by = "Code")

# Normalising
residue_count_df <- residue_count_df %>% mutate(norm_stable_count = (stabilising_count / total_count) * 100,
                                                norm_destable_count = (destabilising_count / total_count) * 100)

# Adding in residue property to colour the bars with 
residue_count_df <- merge(residue_count_df, amino_acids, by = "Code")

# Saving residue counts
write.csv(residue_count_df, file = file.path("Output", "thermodynamics", "residue_analysis", "residue_counts.csv"), row.names = FALSE)

### Stable Plotting ###

# Ordering residue counts from highest to lowest
stable_residue_counts$Code <- factor(x = stable_residue_counts$Code, 
                              levels = stable_residue_counts$Code[order(stable_residue_counts$stabilising_count, decreasing = TRUE)])

residue_count_df$Code <- factor(x = residue_count_df$Code, 
                              levels = residue_count_df$Code[order(residue_count_df$norm_stable_count, decreasing = TRUE)])

ggplot(data = residue_count_df %>% filter(!is.na(norm_stable_count)), 
       aes(x = Code, y = norm_stable_count, fill = property)) +
  geom_col() +
  ylab("Stabilising Occurences / Total Occurences (%)") +
  xlab("") +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Aromatic" = "salmon4", "Hydrophobic" = "goldenrod1", 
                               "Negative" = "dodgerblue", "Polar" = "yellowgreen", 
                               "Positive" = "red3")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        plot.title = element_text(size = 18, face = "bold", colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 5)),
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 20, face = "bold", colour = "black", angle = 90, vjust = 0.5),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 14, colour = "black")) +
  labs(fill = "Property")

ggsave(plot = last_plot(), 
       file = file.path("Output", "thermodynamics", "residue_analysis", "stabilising_residue_type_normalised.png"),
       height = 8, width = 7, bg = "white")

### Destable Plotting ###

# Ordering residue counts from highest to lowest
destable_residue_counts$Code <- factor(x = destable_residue_counts$Code, 
                                     levels = destable_residue_counts$Code[order(destable_residue_counts$destabilising_count, decreasing = TRUE)])

residue_count_df$Code <- factor(x = residue_count_df$Code, 
                                levels = residue_count_df$Code[order(residue_count_df$norm_destable_count, decreasing = TRUE)])

ggplot(data = residue_count_df %>% filter(!is.na(norm_destable_count)), 
       aes(x = Code, y = norm_destable_count, fill = property)) +
  geom_col() +
  ylab("Destabilising Occurences / Total Occurences (%)") +
  xlab("") +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Aromatic" = "salmon4", "Hydrophobic" = "goldenrod1", 
                               "Negative" = "dodgerblue", "Polar" = "yellowgreen", 
                               "Positive" = "red3")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        plot.title = element_text(size = 18, face = "bold", colour = "black", hjust = 0.5),
        axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 5)),
        axis.text.y = element_text(size = 18, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 20, face = "bold", colour = "black", angle = 90, vjust = 0.5),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 14, colour = "black")) +
  labs(fill = "Property")

ggsave(plot = last_plot(), 
       file = file.path("Output", "thermodynamics", "residue_analysis", "destabilising_residue_type_normalised.png"),
       height = 8, width = 7, bg = "white")

# Plotting residue type on line graph

tmp <- filtered_df %>%
  group_by(Pos, Code) %>%
  summarise(mean_energy = mean(total, na.rm=TRUE)) %>%
  ungroup()

tmp2 <- merge(tmp, amino_acids, by = "Code")

# Line Plot
ggplot(data = tmp2, aes(x=Pos, y = mean_energy)) +
  geom_line(linewidth = 0.75) +
  geom_point(aes(colour = property), size = 1.5) +
  #ylab("Mean \u0394G per residue (kcal.mol-1)") +
  ylab(bquote(bold("Mean " * Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
  xlab("Residue Number") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        plot.title = element_text(size = 14, colour = "black", face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"),
        legend.text = element_text(size = 10, colour = "black")) +
  labs(colour = "Property")


# Saving Line plot
ggsave(plot = last_plot(), file = file.path("Output", "thermodynamics", "residue_analysis", "line_plot_with_residue_property.png"), 
       height = 4, width = 10)


############################################################################################################
############################################################################################################
############################################################################################################

####################################
### Sum deltaG per PDB analysis ###
####################################

sum_per_pdb <- filtered_df %>% 
  group_by(pdb_id, Pos) %>%    
  summarise(mean_per_Pos = mean(total, na.rm = TRUE)) %>%
  group_by(pdb_id) %>%
  summarise(total_energy = sum(mean_per_Pos, na.rm = TRUE))

# Loading RMSD Cluster Data
cluster_groups <- read.csv(file.path("Output", "RMSD", "data", "RMSD_cluster_groups.csv"))

# Selecting columns of interest
cluster_groups <- cluster_groups %>%
  select(pdb_id, group)

# Setting group to character
cluster_groups$group <- as.character(cluster_groups$group)

# Remove all instances of duplicated rows based on the pdb_id column
cluster_groups <- cluster_groups %>%
  group_by(pdb_id) %>%
  filter(n() == 1) %>%
  ungroup()


# Merging mean and cluster dataframes
rmsd_cluster_df <- merge(sum_per_pdb, cluster_groups, by = "pdb_id")

# Ordering by group
rmsd_cluster_df$group <- as.numeric(rmsd_cluster_df$group)
rmsd_cluster_df <- rmsd_cluster_df %>% arrange(group)

write.csv(rmsd_cluster_df, file = file.path("Output", "thermodynamics", "deltaG_comparisons", "sum_deltaG_per_PDB_by_RMSD_cluster.csv"),
          row.names = FALSE)

# Plotting
ggplot(data = rmsd_cluster_df, aes(x = factor(group, levels = sort(unique(as.numeric(group)))), y = total_energy)) +
  geom_boxplot(linewidth = 0.75, colour = "black") +
  geom_point(size = 2.5) +
  ylab(bquote(bold(Delta * "G"*degree * " per PDB (kcal.mol"^-1 * ")"))) +
  xlab("RMSD Cluster Group") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 18, face = "bold")) 

ggsave(plot = last_plot(), 
       file = file.path("Output", "thermodynamics", "deltaG_comparisons", "sum_deltaG_per_PDB_by_RMSD_cluster.png"),
       height = 6, width = 8, bg = "white")


### Comparing by fibril source ###

# Loading metadata
selected_pdbs_metadata <- read_delim(file.path("Output", "selected_pdbs_metadata.txt"), 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

# Changing column names
names(selected_pdbs_metadata)[names(selected_pdbs_metadata) == "PDB ID"] <- "PDB"

# creating condition column based on text in fibril origins
selected_pdbs_metadata <- selected_pdbs_metadata %>%
  mutate(condition = case_when(
    str_detect(`Fibril Origins`, regex("extracted|patient|case|atrophy", ignore_case=TRUE)) &
      str_detect(`Fibril Origins`, regex("seed|seeded", ignore_case=TRUE)) ~ "Seeded From Ex Vivo",
    
    str_detect(`Fibril Origins`, regex("extracted|patient|case|atrophy", ignore_case=TRUE)) ~ "Ex Vivo",
    
    !str_detect(`Fibril Origins`, regex("extracted|patient|case|atrophy", ignore_case=TRUE)) &
      str_detect(`Fibril Origins`, regex("seed|seeded", ignore_case=TRUE)) ~ "Seeded From In Vitro",
    
    !str_detect(`Fibril Origins`, regex("extracted|patient|case|atrophy", ignore_case=TRUE)) ~ "In Vitro",
    TRUE ~ "other"
  ))

# Selecting columns of interest
metadata <- selected_pdbs_metadata %>% select(c("PDB", "condition"))


# Creating PDB column
sum_per_pdb <- sum_per_pdb %>%
                  mutate(PDB = if_else(str_detect(pdb_id, "_"), str_extract(pdb_id, "^[^_]+"), pdb_id))

# Merging means with metadata
sum_energy_by_condition <- merge(sum_per_pdb, metadata, by = "PDB")

write.csv(sum_energy_by_condition, file = file.path("Output", "thermodynamics", "deltaG_comparisons", "sum_deltaG_per_PDB_by_fibril_condition.csv"),
          row.names = FALSE)

# Plotting
ggplot(data = sum_energy_by_condition, aes(x = condition, y = total_energy)) +
  geom_boxplot(linewidth = 0.75, colour = "black") +
  geom_point(size = 2.5) +
  #ylab("Mean \u0394G per PDB (kcal.mol-1)") +
  ylab(bquote(bold(Delta * "G"*degree * " per PDB (kcal.mol"^-1 * ")"))) +
  xlab("Fibril Condition") +
  scale_x_discrete(labels = c("Ex vivo" = "Ex Vivo", 
                              "In Vitro" = "In Vitro", 
                              "Seeded From In Vitro" = "Seeded From\nIn Vitro", 
                              "Seeded From Ex Vivo" = "Seeded From\nEx Vivo")) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 18, face = "bold")) 

# Saving plot
ggsave(plot = last_plot(), 
       file = file.path("Output", "thermodynamics", "deltaG_comparisons", "sum_deltaG_per_PDB_by_fibril_condition.png"),
       height = 6, width = 8, bg = "white")


###########################
### Stability Bar Chart ###
###########################

tmp <- filtered_df %>% 
        group_by(Pos) %>%
          summarize(mean_energy = mean(total, na.rm = TRUE)) %>% 
            ungroup() %>%
              mutate(colour_catagory = case_when(
                mean_energy > 0 ~ "above_0",
                mean_energy <= 0 & mean_energy > stabilising_threshold ~ "below_0_above_threshold",
                mean_energy <= stabilising_threshold ~ "below_threshold"
              ))

ggplot(data = tmp, aes(x = Pos, y = mean_energy, fill = colour_catagory)) +
  geom_col() +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = stabilising_threshold, colour = "#f27272", linewidth = 1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0,10000, 10), expand = c(0,0)) + 
  scale_fill_manual(values = c(
    "above_0" = "#6BAF5F",
    "below_0_above_threshold" = "#fabb1b",
    "below_threshold" = "#f27272"
  )) +
  labs(#y = "Mean \u0394G per Residue (kcal.mol-1)",
       y = bquote(bold("Mean " * Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")")),
       x = "Residue") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = "black", face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) 


# Saving plot
ggsave(plot = last_plot(), file = file.path("Output", "thermodynamics", "deltaG_analysis", "bar_chart.png"),
       height = 4, width = 8)
