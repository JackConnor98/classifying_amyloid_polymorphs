cat("\nRunning: thermodynamics_plotting.R\n")

# Installing and Loading Packages
required_packages <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "zoo", "ggrepel", "dendextend", "FSA")
#required_packages <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "zoo", "ggrepel", "dendextend", "FSA", "plotly", "htmlwidgets")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]   

# Install missing packages
if(length(not_installed)) install.packages(not_installed)    

# Not loading plotly and htmlwidgets as they interfere with ggplot's ggsave()
packages_to_load <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "zoo", "ggrepel", "dendextend")

# Load all required packages
invisible(lapply(packages_to_load, require, character.only = TRUE))

# Cleaning Environment
rm(list = ls())


################################################################################

# Creating Directories
if(!dir.exists("Output/thermodynamics/deltaG_analysis")) dir.create("Output/thermodynamics/deltaG_analysis")
if(!dir.exists("Output/thermodynamics/single_PDB_plots")) dir.create("Output/thermodynamics/single_PDB_plots")
if(!dir.exists("Output/thermodynamics/residue_analysis")) dir.create("Output/thermodynamics/residue_analysis")
if(!dir.exists("Output/thermodynamics/pearson_correlation")) dir.create("Output/thermodynamics/pearson_correlation")
if(!dir.exists("Output/thermodynamics/deltaG_comparisons")) dir.create("Output/thermodynamics/deltaG_comparisons")
if(!dir.exists("Output/thermodynamics/stability_clustering")) dir.create("Output/thermodynamics/stability_clustering")

################################################################################

### Importing Data ###

df <- read_csv("Output/thermodynamics/foldx_stability.csv")

# Removing non-cryoEM structures
metadata <- read_delim("Output/selected_pdbs_metadata.txt", 
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
chain_fibril <- read_csv("Output/PDBs/fibrils_extended/chain_fibril.csv")

# Renaming columns to match
names(df)[names(df) == "Mol"] <- "chain"

# Creating new fbril column in df based on matches in chain_fibril
df <- left_join(df, chain_fibril, by = c("PDB", "chain"))

################################################################################

### Removing Head and Tail Chains ###

# Importing data showing which chains are head and tails for each fibril
exterior_chains <- read_csv("Output/PDBs/fibrils_extended/exterior_chains.csv")

filtered_df <- df[!(paste(df$PDB, df$chain) %in% paste(exterior_chains$PDB, exterior_chains$chain)), ]

################################################################################

### Adding in pdb_id to handle fibrils with 2+ distinct polymorphs ###

# Importing data
COM_and_fibril <- read.csv("Output/PDBs/COM_and_fibril.csv")

# Selecting columns of interest
pdb_id_data <- COM_and_fibril %>% select(c("PDB", "pdb_id", "fibril"))

# Removing duplicate rows
pdb_id_data <- pdb_id_data %>% distinct()

# Adding pdb_id to filtered_df
filtered_df <- merge(pdb_id_data, filtered_df, by = c("PDB", "fibril"))

#################################################################################################################

######################################################################
### Analysing all data including single residues with low Q-scores ###
######################################################################

# Calculating the mean delta G per residue for each PDB 
no_q_score_filter <- filtered_df %>%
                      group_by(pdb_id, Pos) %>%
                        summarize(mean_energy = mean(total, na.rm = TRUE)) %>%
                          ungroup()
 
# Sliding window average (window size 5)
no_q_score_filter <- no_q_score_filter %>%
            group_by(pdb_id) %>%
              mutate(sliding_window = rollmean(mean_energy, k = 5, FUN = mean, align = "center", fill = NA)) %>%
                ungroup()

# Scaling data so mean = 0 and SD = 1
no_q_score_filter <- no_q_score_filter %>%
                      group_by(pdb_id) %>%
                        mutate(scaled = scale(sliding_window)[,1]) %>%
                          ungroup()

write.csv(no_q_score_filter, file = "Output/thermodynamics/foldx_stability_unfiltered_normalised.csv", row.names = FALSE)

# Plotting unfiltered data
ggplot(data = no_q_score_filter, aes(x=Pos)) +
  geom_line(aes(y = scaled, group = interaction(pdb_id, cumsum(c(0, diff(Pos) != 1)))), 
            linewidth = 0.75, alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 0.75, linetype = "solid", alpha = 0.9) +
  #ylab("Mean Scaled \u0394G per residue (kcal.mol-1)") +
  ylab(bquote(bold("Mean Scaled " * Delta * "G" * degree * " per residue (kcal.mol"^-1 * ")"))) +
  xlab("Residue Number") +
  scale_x_continuous(breaks = seq(0, 99999, by = 1),
                    labels = ifelse(seq(0, 99999, by = 1) %% 10 == 0, seq(0, 99999, by = 1), ""),
                    expand = expansion(mult = c(0.01,0.01))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        #panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text = element_text(size = 14, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 1.5), 
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.position = "none") 

ggsave(plot = last_plot(), file = "Output/thermodynamics/deltaG_analysis/stability_line_plot_sliding_window_unfiltered_scaled.png",
       height = 6, width = 12)


################################################################################

##############################################
### Removing residues with poor resolution ###
##############################################

# Importing data
high_resolution_residues <- read.csv("Output/Validation/high_resolution_residues.csv")

# Renaming columns to match df
names(high_resolution_residues)[names(high_resolution_residues) == "resno"] <- "Pos"

# Merging with main data frame
filtered_df <- merge(filtered_df, high_resolution_residues, by = c("pdb_id", "fibril", "Pos"))

##############################################################
### Removing only low residue PDBs not individual residues ###
##############################################################

# Getting a list of PDBs with high resolution
#high_resolution_PDBs <- unique(high_resolution_residues$pdb_id)

# Keeping only PDBs found in high_resolution_PDBs
#filtered_df <- filtered_df %>% filter(pdb_id %in% high_resolution_PDBs)

################################################################################

# Saving filtered_df
write.csv(filtered_df, file = "Output/thermodynamics/foldx_stability_filtered.csv", row.names = FALSE)

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
                    summarize(mean_energy = mean(total, na.rm = TRUE))

# Saving Data
#write.csv(mean_per_res, file = "Output/thermodynamics/mean_deltaG_per_residue.csv", row.names = FALSE)


# Calculating the threshold for residues to be considered stabilising
mean_stability <- mean(mean_per_res$mean_energy, na.rm = TRUE)
SD_stability <- sd(mean_per_res$mean_energy, na.rm = TRUE)
stabilising_threshold <- mean_stability - SD_stability
destabilising_threshold <- mean_stability + SD_stability


# Plotting

### Density ###
ggplot() +
  geom_density(data = mean_per_res, aes(x = mean_energy, group = pdb_id), 
               colour=alpha("black", 0.25), linewidth = 1) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 1, linetype = "dashed", alpha = 0.9) +
  ylab("Density") +
  #xlab("\u0394G per residue (kcal.mol-1)") +
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

ggsave(plot = last_plot(), file = "Output/thermodynamics/deltaG_analysis/stability_density_plot.png",
       height = 6, width = 8)


### Annotating residues with DeltaG <= stabilising_threshold ###


# Function to find consecutive ranges
find_ranges <- function(x) {
  if(length(x) == 0) return(character(0))
  x <- sort(x)
  ranges <- c()
  start <- x[1]
  end <- x[1]

  for(i in 2:length(x)) {
    if(x[i] == end + 1) {
      end <- x[i]
    } else {
      if(start == end) {
        ranges <- c(ranges, as.character(start))
      } else {
        ranges <- c(ranges, paste0(start, "-", end))
      }
      start <- x[i]
      end <- x[i]
    }
  }

  if(start == end) {
    ranges <- c(ranges, as.character(start))
  } else {
    ranges <- c(ranges, paste0(start, "-", end))
  }

  return(ranges)
}


# Ungroup the data first
mean_per_res <- mean_per_res %>% ungroup()

# Filter the points where mean_energy is < stabilising_threshold
annotation_points <- mean_per_res %>%
  filter(mean_energy <= stabilising_threshold) %>%
  select(Pos, mean_energy) %>%  
  distinct(Pos, .keep_all = TRUE)


# Generate labels with grouped consecutive ranges
grouped_labels <- find_ranges(annotation_points$Pos)

# Create a dataframe for labels
labels_df <- data.frame(
  Pos = sapply(grouped_labels, function(x) {
    if(grepl("-", x)) {
      return(mean(as.numeric(unlist(strsplit(x, "-")))))
    } else {
      return(as.numeric(x))
    }
  }),
  Label = grouped_labels,
  mean_energy = sapply(grouped_labels, function(x) {
    if(grepl("-", x)) {
      range_values <- as.numeric(unlist(strsplit(x, "-")))
      return(mean(annotation_points$mean_energy[annotation_points$Pos %in% range_values]))
    } else {
      return(annotation_points$mean_energy[annotation_points$Pos == as.numeric(x)])
    }
  })
)

# Removing asyn incorrectly modelled residue
#tmp <- mean_per_res %>% filter(pdb_id != "7v4a" | Pos != 35)

p <- ggplot(data = mean_per_res, aes(x=Pos)) +
  geom_line(aes(y = mean_energy, group = interaction(pdb_id, cumsum(c(0, diff(Pos) != 1)))), 
            linewidth = 0.75, alpha = 0.25) +
  geom_hline(yintercept = 0, colour = "blue", linewidth = 0.75, linetype = "solid", alpha = 0.9) +
  
  #geom_hline(yintercept = destabilising_threshold, colour = "blue", linewidth = 1, linetype = "solid", alpha = 0.9) +
  #geom_hline(yintercept = stabilising_threshold, colour = "blue", linewidth = 1, linetype = "solid", alpha = 0.9) +
  
  #geom_label_repel(data = labels_df,
  #                 aes(x = Pos, y = mean_energy, label = Label, fill = after_scale(alpha("grey", 0.3))),
  #                 colour = "black",
  #                 box.padding  = 0.5,
  #                 point.padding = 0.2,
  #                 segment.color = "black",
  #                 segment.size = 1.5,
  #                 label.size = 0.8,
  #                 size = 6,
  #                 max.overlaps = 100) +
  
  #ylab("\u0394G per residue (kcal.mol-1)") +
  ylab(bquote(bold(Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
  xlab("Residue Number") +
  scale_x_continuous(breaks = seq(0, 99999, by = 1), 
                     labels = ifelse(seq(0, 99999, by = 1) %% 10 == 0, seq(0, 99999, by = 1), ""),
                     expand = expansion(mult = c(0.01,0.01))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        #panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text = element_text(size = 14, colour = "black", face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold", vjust = 1.5), 
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.position = "none") #+ ylim(-3,5)

plot(p)

ggsave(plot = p, file = "Output/thermodynamics/deltaG_analysis/stability_line_plot.png",
       height = 6, width = 12)

# Convert ggplot to interactive plotly
#interactive_plot <- plotly::ggplotly(p)

# Save interactive plot as a HTML file
#htmlwidgets::saveWidget(interactive_plot, file = "Output/thermodynamics/deltaG_analysis/stability_line_plot_interactive.html")


####################################
####################################
##### Sliding window averaging #####
####################################
####################################

window_size <- 5

sliding_window <- mean_per_res %>%
  group_by(pdb_id) %>%
  mutate(sliding_window = rollapply(mean_energy, width = window_size, FUN = mean, align = "center", fill = NA))

# Saving Data
write.csv(sliding_window, file = "Output/thermodynamics/foldx_stability_filtered_normalised.csv", row.names = FALSE)


# Ungroup the data first
sliding_window <- sliding_window %>% ungroup()


# Filter the points where sliding_window is < -0.5
annotation_points <- sliding_window %>%
                      filter(sliding_window <= -0.5) %>%
                        select(Pos, sliding_window) %>%  
                          distinct(Pos, .keep_all = TRUE)



# Generate labels with grouped consecutive ranges
grouped_labels <- find_ranges(annotation_points$Pos)

# Create a dataframe for labels
labels_df <- data.frame(
  Pos = sapply(grouped_labels, function(x) {
    if(grepl("-", x)) {
      return(mean(as.numeric(unlist(strsplit(x, "-")))))
    } else {
      return(as.numeric(x))
    }
  }),
  Label = grouped_labels,
  sliding_window = sapply(grouped_labels, function(x) {
    if(grepl("-", x)) {
      range_values <- as.numeric(unlist(strsplit(x, "-")))
      return(mean(annotation_points$sliding_window[annotation_points$Pos %in% range_values]))
    } else {
      return(annotation_points$sliding_window[annotation_points$Pos == as.numeric(x)])
    }
  })
)
 
p <- ggplot(data = sliding_window, aes(x=Pos)) +
 geom_line(aes(y = sliding_window, group = interaction(pdb_id, cumsum(c(0, diff(Pos) != 1)))), 
           linewidth = 1, alpha = 0.25) +
 geom_hline(yintercept = 0, colour = "blue", linewidth = 1, linetype = "solid", alpha = 0.9) +
 #geom_label_repel(data = labels_df,
#                  aes(x = Pos, y = sliding_window, label = Label, fill = after_scale(alpha("grey", 0.3))),
#                  colour = "black",
#                  box.padding  = 0.5, 
#                  point.padding = 0.2,
#                  segment.color = "black",
#                  segment.size = 1.5,
#                  label.size = 0.8,
#                  size = 6) +
 ggtitle("Sliding Window of 5 Mean") +
 #ylab("\u0394G per residue (kcal.mol-1)") +
 ylab(bquote(bold(Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
 xlab("Residue Number") +
 scale_x_continuous(breaks = seq(0,99999, by = 10)) +
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
       legend.position = "none") 

plot(p)

ggsave(plot = p, file = "Output/thermodynamics/deltaG_analysis/stability_line_plot_sliding_window.png",
      height = 6, width = 12)

# Convert ggplot to interactive plotly
#interactive_plot <- plotly::ggplotly(p)

# Save interactive plot as a HTML file
#htmlwidgets::saveWidget(interactive_plot, file = "Output/thermodynamics/deltaG_analysis/stability_line_plot_sliding_window_interactive.html")


############################################################################################

### Plotting each PDB individually so I can make a movie and visualize the shape of each ###

############################################################################################

pdb_names <- unique(sliding_window$pdb_id)

for (pdb in pdb_names) {

  x <- sliding_window %>% filter(pdb_id %in% pdb)

  p <- ggplot(data = x, aes(x=Pos)) +
    geom_line(aes(y=sliding_window, group = cumsum(c(0, diff(Pos) != 1))), linewidth = 2) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 1, linetype = "dashed", alpha = 0.9) +
    geom_hline(yintercept = -0.5, colour = "red", linewidth = 1, linetype = "dashed", alpha = 0.9) +
    ggtitle(paste(pdb, "Sliding Window of 5 Mean", sep = " - ")) +
    #ylab("\u0394G per residue (kcal.mol-1)") +
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
  
  if ("8azs" %in% pdb_names) {
    p + scale_x_continuous(limits = c(0, 42), breaks = seq(0,99999, by = 10))
  } 
  if ("6cu7" %in% pdb_names) {
    p + scale_x_continuous(limits = c(0, 140), breaks = seq(0,99999, by = 10)) 
  }
  
  p

  ggsave(plot = last_plot(), file = paste("Output/thermodynamics/single_PDB_plots/", pdb, "_stability.png"),
         height = 6, width = 12)

}

################################################################################

### Clustering by per residue stability ### 

################################################################################

wide_df <- mean_per_res %>% pivot_wider(names_from = "Pos",
                                        values_from = "mean_energy")

# Removing reference columns
cluster_data <- subset(wide_df, select = setdiff(names(wide_df), c("pdb_id")))

# Replacing NAs with Mean + 3SD (Occurs when PDBs have no overlapping residues)
cluster_data[is.na(cluster_data)] <- mean(as.matrix(cluster_data), na.rm = TRUE) + 3*sd(as.matrix(cluster_data), na.rm = TRUE)

# Scaling Cluster Data (makes mean of all columns = 0 and SD = 0-1)
scaled_cluster_data <- as.data.frame(scale(cluster_data))
summary(scaled_cluster_data)

# Calculating Distance
cluster_dist <- dist(scaled_cluster_data, method = "euclidean")

# Clustering
hc_average <- hclust(cluster_dist, method = "average")

# Making labels
pdb_id_labels <- wide_df$pdb_id[hc_average$order]

# Setting cut_height to the mean distance
cut_height <- mean(cluster_dist)

# Plotting
avg_dend_obj <- as.dendrogram(hc_average)
avg_col_dend <- color_branches(avg_dend_obj, h = cut_height)
labels(avg_col_dend) <- pdb_id_labels
plot(avg_col_dend, ylab = "Euclidean Distance") # horiz = TRUE
abline(h = cut_height, col = "red", lty = 2, lw = 2)


# Save the cluster plot as a PNG file
png("Output/thermodynamics/stability_clustering/foldx_stability_cluster_dendrogram.png", width = 1600, height = 500, res = 100)
plot(avg_col_dend, ylab = "Euclidean Distance", xlab = "PDB ID")
abline(h = cut_height, col = "red", lty = 2, lw = 2)
title("Cluster Dendrogram - Per Residue FoldX Stability")
dev.off()

# Extracting pdb_id from each cluster group
clusters <- cutree(hc_average, h = cut_height)
grouped_pdbs <- split(wide_df$pdb_id, clusters)

# Convert the list to a dataframe
cluster_groups <- do.call(rbind, lapply(names(grouped_pdbs), function(name) {
                    data.frame(group = name, pdb_id = grouped_pdbs[[name]], stringsAsFactors = FALSE)
                  }))
# Saving Cluster Groups
write.csv(cluster_groups, file = "Output/thermodynamics/stability_clustering/stability_cluster_groups.csv", row.names = FALSE)



############################################################################################

### Calculating the Pearson Correlation Coefficient Between PDB Stability ### 

############################################################################################ 

#tmp <- read_csv("Output/thermodynamics/foldx_stability_unfiltered_normalised.csv")

wide_df <- sliding_window %>%
  select(c("pdb_id", "Pos", "sliding_window")) %>%
  pivot_wider(names_from = pdb_id, values_from = sliding_window) 

# Select columns for correlation calculation (excluding the first column)
columns_to_correlate <- wide_df[, -1]

# Calculate correlation matrix, ignoring NA values
correlation_matrix <- cor(columns_to_correlate, use = "pairwise.complete.obs")

# Convert correlation matrix to data frame
correlation_df <- as.data.frame(as.table(correlation_matrix))
names(correlation_df) <- c("Var1", "Var2", "Correlation")

# Ordering heatmap by clustering
correlation_df$Var1 <- factor(correlation_df$Var1, levels = pdb_id_labels)
correlation_df$Var2 <- factor(correlation_df$Var2, levels = pdb_id_labels)

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

ggsave(plot = last_plot(), file = paste("Output/thermodynamics/pearson_correlation/Pearson_heatmap.png"),
       height = 11, width = 12, bg = "white")

write.csv(correlation_df, file = "Output/thermodynamics/pearson_correlation/Pearson_Scores.csv", row.names = FALSE)


### Correlation Box Plot ###

# Adding a dummy column to use for x-axis
correlation_df$dummy <- ""

# Creating a comparison column
correlation_df <- correlation_df %>%
                    mutate(comparison = paste0(Var1, "-", Var2))

# Removing Self Comparisons
correlation_df <- correlation_df %>% filter(Var1 != Var2)

# Define a function to identify outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm = TRUE) - 2 * IQR(x, na.rm = TRUE) | 
         x > quantile(x, 0.75, na.rm = TRUE) + 2 * IQR(x, na.rm = TRUE))
}

# Apply the function to your data frame
correlation_df <- correlation_df %>%
  mutate(outlier = ifelse(is_outlier(Correlation), as.character(comparison), NA))

# plot
ggplot(data = correlation_df, aes(x = dummy, y = Correlation)) +
  geom_violin(fill = "grey70", linewidth = 1, width = 0.5) +
  geom_boxplot(size = 1, colour = "black", width = 0.1, fill = NA) +
  #geom_text_repel(aes(label = outlier), size = 3, fontface = "bold", 
  #                na.rm = TRUE, hjust = -0.15, vjust = 0.5) + 
  scale_y_continuous(name = "Pearson Correlation Between\nPDBs Per Residue FoldX Stability") +
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


ggsave(plot = last_plot(), file = paste("Output/thermodynamics/pearson_correlation/Pearson_boxplot.png"),
       height = 8, width = 8)

### Saving Median and Mean Correlation ###
sink("Output/thermodynamics/pearson_correlation/average_correlation.txt")
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


# Sliding window averaging 

window_size <- 5

sliding_window <- pdb_means %>%
  group_by(pdb_id) %>%
  mutate(across(starts_with("mean_"), 
                ~ rollapply(.x, width = window_size, FUN = mean, align = "center", fill = NA))) %>%
  ungroup()


# Reshape the dataframe from wide to long format
#long_df <- sliding_window %>%
#  pivot_longer(cols = starts_with("mean_"), names_to = "variable", values_to = "value") 

long_df <- pdb_means %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "variable", values_to = "value") 


# Removing NA values
long_df <- na.omit(long_df)

# Get the order of columns in sliding_window
column_order <- names(sliding_window)[-1]  # Exclude the Pos column

# Remove "mean_" pattern
long_df$variable <- gsub("mean_", "", long_df$variable)
column_order <- gsub("mean_", "", column_order)

# Setting order of long_df (important for cumsum() to calculate line breaks)
long_df <- long_df[order(long_df$variable, long_df$pdb_id, long_df$Pos), ]

# Setting plot order
long_df$variable <- factor(long_df$variable, levels = column_order)

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

ggsave(plot = last_plot(), file = "Output/thermodynamics/deltaG_analysis/all_energies_separate_pdb.png",
       height = 12, width = 15, bg = "white")

### Averaged per Residue ###

#selected_columns <- c("total", "backHbond", "sideHbond", "energy_VdW", "electro", "energy_SolvP", "energy_SolvH", "energy_vdwclash", "entrop_sc", "entrop_mc")

selected_columns <- colnames(filtered_df[11:32])

# Calculate the mean of the selected columns and store it in new columns
means <- filtered_df %>%
  group_by(Pos) %>%
  summarise(across(all_of(selected_columns), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"))

# Sliding window averaging 
window_size <- 5

sliding_window <- means %>%
  mutate(across(starts_with("mean_"), 
                ~ rollapply(.x, width = window_size, FUN = mean, align = "center", fill = NA))) %>%
  ungroup()

# Reshape the dataframe from wide to long format
long_df <- sliding_window %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "variable", values_to = "value")

# Removing NA values
long_df <- na.omit(long_df)

# Get the order of columns in sliding_window
column_order <- names(sliding_window)[-1]  # Exclude the Pos column

# Remove "mean_" pattern
long_df$variable <- gsub("mean_", "", long_df$variable)
column_order <- gsub("mean_", "", column_order)

# Setting order of long_df$variable
long_df$variable <- factor(x=long_df$variable, levels = column_order)

# Setting plot order
long_df$variable <- factor(long_df$variable, levels = column_order)

# Create the plot
ggplot(long_df, aes(x = Pos, y = value, color = variable, group = variable)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  facet_wrap(~variable) +
  #labs(x = "Position", y = "\u0394G per Residue (kcal.mol-1)", title = "Mean Scores - All PDBs") +
  
  labs(title = "Mean Scores - All PDBs", x = "Position",
       y = bquote(bold("Mean Scaled " * Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
  
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.line = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

ggsave(plot = last_plot(), file = "Output/thermodynamics/deltaG_analysis/all_energies_mean_pdb.png",
       height = 12, width = 15, bg = "white")


#####################################################################################################
#####################################################################################################
#####################################################################################################

### Characterizing Residue Position ### 

#####################################################################################################
#####################################################################################################
#####################################################################################################

filtered_df$Code[filtered_df$Code == "H1S"] <- "HIS"
filtered_df$Code[filtered_df$Code == "H2S"] <- "HIS"

# Counting how many PDBs each residue position is found in
all_pos_counts <- filtered_df %>% count(Pos)

# Counting how many PDBs each residue position is found to be stabilising
stabilising_df <- filtered_df %>% filter(total <= stabilising_threshold)
stable_pos_counts <- stabilising_df %>% count(Pos)

# Counting how many PDBs each residue position is found to be destabilising
destabilising_df <- filtered_df %>% filter(total >= destabilising_threshold)
destable_pos_counts <- destabilising_df %>% count(Pos)



# Normalising the stablising counts by the total counts
names(all_pos_counts)[names(all_pos_counts) == "n"] <- "total_count"
names(stable_pos_counts)[names(stable_pos_counts) == "n"] <- "stabilising_count"
names(destable_pos_counts)[names(destable_pos_counts) == "n"] <- "destabilising_count"

# Merging Data
pos_count_df <- merge(stable_pos_counts, destable_pos_counts, by = "Pos", all = TRUE)
pos_count_df <- merge(pos_count_df, all_pos_counts, by = "Pos")

# Normalising Counts
pos_count_df <- pos_count_df %>% mutate(norm_stable_count = (stabilising_count / total_count) * 100,
                                        norm_destable_count = (destabilising_count / total_count) * 100)


### Stable Plotting ###

ggplot() +
  geom_col(data = pos_count_df, aes(x = Pos, y = norm_stable_count)) +
  ggtitle("Normalised by Total Number of Occurences") +
  ylab("Number of Stabilising Occurences") +
  xlab("Residue Number") +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,9999,10)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold", colour = "black", hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", colour = "black"))

ggsave(plot = last_plot(), 
       file = "Output/thermodynamics/residue_analysis/stabilising_residue_number_normalised.png",
       height = 6, width = 8, bg = "white")

### Destable Plotting ###

ggplot() +
  geom_col(data = pos_count_df, aes(x = Pos, y = norm_destable_count)) +
  ggtitle("Normalised by Total Number of Occurences") +
  ylab("Number of Destabilising Occurences") +
  xlab("Residue Number") +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,9999,10)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold", colour = "black", hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", colour = "black"))

ggsave(plot = last_plot(), 
       file = "Output/thermodynamics/residue_analysis/destabilising_residue_number_normalised.png",
       height = 6, width = 8, bg = "white")


#####################################################################################################
#####################################################################################################
#####################################################################################################

### Characterizing Residue Type ### 

#####################################################################################################
#####################################################################################################
#####################################################################################################

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
       file = "Output/thermodynamics/residue_analysis/stabilising_residue_type_normalised.png",
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
       file = "Output/thermodynamics/residue_analysis/destabilising_residue_type_normalised.png",
       height = 8, width = 7, bg = "white")

############################################################################################################

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
ggsave(plot = last_plot(), file = "Output/thermodynamics/residue_analysis/line_plot_with_residue_property.png", height = 4, width = 10)


############################################################################################################
############################################################################################################
############################################################################################################

####################################
### Sum deltaG per PDB analysis ###
####################################

# Reimporting FSA library as it doesn't work when imported at the start for some reason
library(FSA)

# Calculating mean deltaG for each PDB
sum_per_pdb <- filtered_df %>%
  group_by(pdb_id) %>%
  #summarize(mean_energy = mean(total, na.rm = TRUE))
  summarize(total_energy = sum(total, na.rm = TRUE)) # Summing is better than mean as it takes length into account



# Loading RMSD Cluster Data
cluster_groups <- read.csv("Output/RMSD/RMSD_cluster_groups.csv")

# Removing _chainID to just use PDB codes 
# cluster_groups <- cluster_groups %>%
#   mutate(pdb_id = sapply(str_split(pdb_id, "_"), function(x) x[1]))

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

write.csv(rmsd_cluster_df, file = "Output/thermodynamics/deltaG_comparisons/sum_deltaG_per_PDB_by_RMSD_cluster.csv",
          row.names = FALSE)

if (length(unique(rmsd_cluster_df$group)) > 1) {
  
  # Kruskal-Wallis Test
  result <- kruskal.test(total_energy ~ group, data = rmsd_cluster_df)
  
  print(result)
  
  # Making group column a factor for dunnTest()
  rmsd_cluster_df$group <- as.factor(rmsd_cluster_df$group)
  
  # Dunn Test (p-value correction)
  dunn_test <- dunnTest(total_energy ~ group, data = rmsd_cluster_df, method = "holm")
  
  # Extract results and reorder based on the numeric values in the Comparison column
  sorted_dunn_test <- dunn_test$res %>%
    mutate(Comparison = gsub(" - ", "_", Comparison)) %>%  # Replace " - " with "_"
    separate(Comparison, into = c("Group1", "Group2"), sep = "_", convert = TRUE) %>%  # Split into two numeric columns
    arrange(Group1, Group2)  # Sort properly
  
  print(sorted_dunn_test)
  
  # Saving Tukey's HSD Results
  sink("Output/thermodynamics/deltaG_comparisons/RMSD_cluster_stats.txt")
  #print(tukey_result)
  print(sorted_dunn_test)
  sink()
  
} else {
  print("Only one condition group - Unable to carry out statistical analysis")
}

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


  # # Adding stats -- ASYN
  # annotate("segment", x = 3.02, xend = 3.98, y = 900, yend = 900, linewidth = 1.5) +
  # annotate("text", x = 3.5, y = 960, label = "NS", size = 6, fontface = "bold") +
  # annotate("segment", x = 4.02, xend = 5.98, y = 1000, yend = 1000, linewidth = 1.5) +
  # annotate("text", x = 5, y = 1060, label = "NS", size = 6, fontface = "bold") +
  # annotate("segment", x = 6.02, xend = 6.98, y = 900, yend = 900, linewidth = 1.5) +
  # annotate("text", x = 6.5, y = 960, label = "NS", size = 6, fontface = "bold") +
  # annotate("segment", x = 3.02, xend = 6.98, y = 1250, yend = 1250, linewidth = 1.5) +
  # annotate("text", x = 5, y = 1310, label = "NS", size = 6, fontface = "bold")


ggsave(plot = last_plot(), 
       file = "Output/thermodynamics/deltaG_comparisons/sum_deltaG_per_PDB_by_RMSD_cluster.png",
       height = 6, width = 8, bg = "white")


### Trying with in vitro vs ex vivo ###

# Loading metadata
selected_pdbs_metadata <- read_delim("Output/selected_pdbs_metadata.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

# Changing column names
names(selected_pdbs_metadata)[names(selected_pdbs_metadata) == "PDB ID"] <- "PDB"

# creating condition column based on text in fibril origins
selected_pdbs_metadata <- selected_pdbs_metadata %>%
  mutate(condition = case_when(
    str_detect(`Fibril Origins`, regex("extracted|patient|case|atrophy", ignore_case=TRUE)) &
      str_detect(`Fibril Origins`, regex("seed|seeded", ignore_case=TRUE)) ~ "Seeded From Ex Vivo",
    
    str_detect(`Fibril Origins`, regex("extracted|patient|case|atrophy", ignore_case=TRUE)) ~ "Ex vivo",
    
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

write.csv(sum_energy_by_condition, file = "Output/thermodynamics/deltaG_comparisons/sum_deltaG_per_PDB_by_fibril_condition.csv",
          row.names = FALSE)

if (length(unique(sum_energy_by_condition$condition)) > 1) {
  
  # Kruskal-Wallis Test
  result <- kruskal.test(total_energy ~ condition, data = sum_energy_by_condition)
  
  print(result)
  
  # Making condition into a factor for dunnTest()
  sum_energy_by_condition$condition <- as.factor(sum_energy_by_condition$condition)
  
  # Dunn Test (p-value correction)
  dunn_test <- dunnTest(total_energy ~ condition, data = sum_energy_by_condition, method = "holm")
  
  # Extract results and reorder based on the numeric values in the Comparison column
  sorted_dunn_test <- dunn_test$res %>%
    mutate(Comparison = gsub(" - ", "_", Comparison)) %>%  # Replace " - " with "_"
    separate(Comparison, into = c("Condition1", "Condition2"), sep = "_", convert = TRUE) %>%  # Split into two numeric columns
    arrange(Condition1, Condition2)  # Sort properly
  
  print(sorted_dunn_test)
  
  
  # Saving Tukey's HSD Results
  sink("Output/thermodynamics/deltaG_comparisons/fibril_condition_stats.txt")
  print(sorted_dunn_test)
  sink()
  
} else {
  print("Only one condition group - Unable to carry out statistical analysis")
}

# Splitting seeded onto two lines to allow larger font sizes


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

  # # Adding stats -- ASYN
  # annotate("segment", x = 1.02, xend = 1.98, y = 1350, yend = 1350, linewidth = 1.5) +
  # annotate("text", x = 1.5, y = 1420, label = "P.adj = 0.001", size = 5) +
  # annotate("segment", x = 1.02, xend = 2.98, y = 1500, yend = 1500, linewidth = 1.5) +
  # annotate("text", x = 2, y = 1570, label = "P.adj = 0.026", size = 5)

  # # Adding stats -- TAU
  # annotate("segment", x = 1.02, xend = 1.98, y = 1480, yend = 1480, linewidth = 1.5) +
  # annotate("text", x = 1.5, y = 1550, label = "P.adj = 0.320", size = 5)


  # Saving plot
ggsave(plot = last_plot(), 
       file = "Output/thermodynamics/deltaG_comparisons/sum_deltaG_per_PDB_by_fibril_condition.png",
       height = 6, width = 8, bg = "white")


##################################################################################################
##################################################################################################
##################################################################################################

########################
### Energy Landscape ###
########################

##################################################################################################
##################################################################################################
##################################################################################################

# Improting RMSD Cluster Dendrogram Order
dendrogram_order <- read_csv("Output/RMSD/dendrogram_order.csv")

# Loading RMSD Cluster Data
cluster_groups <- read.csv("Output/RMSD/RMSD_cluster_groups.csv")

# Setting group to character
cluster_groups$group <- as.character(cluster_groups$group)

# Merging dataframes
tmp <- merge(merge(sum_energy_by_condition, dendrogram_order, by = "pdb_id"),
                                             cluster_groups, by = "pdb_id")



ggplot(data = tmp, aes(x = order, y = total_energy)) +
  geom_point(aes(colour = group, shape = condition), size = 4, stroke = 1.5) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "black", linewidth = 1.5) +
  scale_colour_brewer(palette = "Set1") +
  scale_shape_manual(values = c(15,16,0,1)) +
  labs(title = "",
       x = "RMSD Cluster Order",
       #y = "Mean \u0394G per PDB (kcal.mol-1)") +
       y = bquote(bold("Mean " * Delta * "G"*degree * " per PDB (kcal.mol"^-1 * ")")))+
  
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text = element_text(size = 14, colour = "black", face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black")) +
  labs(colour = "Cluster Group",
       shape = "Fibril Origins")

ggsave(plot = last_plot(), file = "Output/thermodynamics/deltaG_comparisons/deltaG_vs_RMSD_cluster_order.png",
       height = 6, width = 16)


##################################################################################################
##################################################################################################
##################################################################################################

###########################
### Stability Bar Chart ###
###########################

##################################################################################################
##################################################################################################
##################################################################################################

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
ggsave(plot = last_plot(), file = "Output/thermodynamics/deltaG_analysis/bar_chart.png",
       height = 4, width = 8)












