cat("Running: RMSD_analysis.R\n")

############################################################################################

# Installing and Loading Packages
required_packages <- c("dplyr", "stringr", "tidyr","readr", "ggplot2", "dendextend", "RColorBrewer", "tibble")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]

# Install missing packages
if(length(not_installed)) install.packages(not_installed)

# Load all required packages
invisible(lapply(required_packages, require, character.only = TRUE))
  
# Cleaning Environment
rm(list = ls())

############################################################################################

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

# Convert 'custom_cut_height' parameter to numeric
custom_cut_height <- as.numeric(config$custom_cut_height)

# Check the custom_cut_height
print(paste("Custom cut height:", custom_cut_height))

############################################################################################

output_dir <- file.path("Output", "RMSD")

data_path <- file.path("Output", "RMSD", "data")

single_pdb_dir <- file.path("Output", "RMSD", "single_reference")

# Creating Directories
if(!dir.exists(single_pdb_dir)) dir.create(single_pdb_dir)

############################################################################################

# Read in df
df <- read.csv(file.path(data_path, "pairwise_rmsd.csv"))

############################################################################################

### Removing Low Resolution Residues ###

# Importing data
high_resolution_residues <- read.csv("Output/Validation/high_resolution_residues.csv")

high_res_pdb_ids <- unique(high_resolution_residues$pdb_id)

filtered_df <- df %>% filter(ref_name %in% high_res_pdb_ids) %>%
  filter((mob_name %in% high_res_pdb_ids))

############################################################################################

# Getting amyloid names from first chain only
amyloid_names <- unique(filtered_df$ref_name)

for (i in amyloid_names) {
    
  x <- filtered_df %>% filter(ref_name == i)
  
  # Plotting RMSD for each reference
  ggplot(x, aes(x = reorder(mob_name, rmsd), y = rmsd)) +
    geom_point(size = 4) +
    ggtitle(paste(i, " - Chain Comparison", sep = "")) +
    ylab("RMSD") +
    xlab("PDB") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(filtered_df$rmsd) * 1.05),
                       breaks = seq(0,1000,5)) +
    theme_bw() +
    theme(panel.grid.major.y = element_line(colour = "grey30", linewidth = 0.25,
                                            linetype = "dashed"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
          axis.title = element_text(size = 18, face = "bold", colour = "black"),
          axis.text.x = element_text(size = 13, colour = "black",
                                     angle = 90, vjust = 0.4, hjust = 1),
          axis.text.y = element_text(size = 14, colour = "black")) 
    
  ggsave(plot = last_plot(), file = file.path(single_pdb_dir, paste0(i, "_RMSD.png")), 
         width = 20, height = 8)
  
}

############################################################################################

# Wide format data
df_wide <- pivot_wider(
  filtered_df,
  names_from = mob_name,
  values_from = rmsd
)

names(df_wide)[names(df_wide) == "ref_name"] <- "pdb_id"

# Removing reference columns
cluster_data <- subset(df_wide, select = setdiff(names(df_wide), c("pdb_id")))

# Setting NA values to the mean + 3SD - This occurs when two PDBs have no overlapping residues 
NA_distance <- mean(as.matrix(cluster_data), na.rm = TRUE) + 3*sd(as.matrix(cluster_data), na.rm = TRUE)

cluster_data[is.na(cluster_data)] <- NA_distance

# Scaling Cluster Data (makes mean of all columns = 0 and SD = 0-1)
scaled_cluster_data <- as.data.frame(scale(cluster_data))

# Saving cluster data
write.csv(scaled_cluster_data, file = file.path(data_path, "scaled_cluster_data.csv"), row.names = FALSE)

############################################################################################

################################
### Clustering Based On RMSD ###
################################

# Calculating Distance
cluster_dist <- dist(scaled_cluster_data, method = "euclidean")

# Clustering
hc_average <- hclust(cluster_dist, method = "average")


### Scree Plot to assigning cut height ###
scree_data <- hc_average$height %>%
              as_tibble() %>%
              add_column(groups = length(hc_average$height):1) %>%
              rename(height=value)

# Saving Scree data
write.csv(scree_data, file = file.path(data_path, "scree_data.csv"), row.names = FALSE)

# Plotting 
ggplot(scree_data, aes(x=groups, y=height)) +
       geom_point(size = 1.5) +
       geom_line(linewidth = 0.5) +
       annotate(geom = "segment", linewidth = 0.5,, colour = "red",
                x = -Inf, xend = Inf,
                y = custom_cut_height, yend = custom_cut_height) +
  scale_x_continuous(expand = c(0,0.75), breaks = seq(0,9999,10)) +
  scale_y_continuous(#limits = c(0,30),
                     expand = c(0,0.75)) +
  labs(x = "Groups", y = "Euclidean Distance") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.5, linetype = "dashed"),
        axis.title = element_text(size = 16, colour = "black", face = "bold"),
        axis.text = element_text(size = 14, colour = "black"))

# Saving the scree plot
ggsave(plot = last_plot(), filename = file.path(output_dir, "Scree_plot.png"), height = 4, width = 8)



# Making labels
pdb_id_labels <- df_wide$pdb_id[hc_average$order]

# Creating a dataframe to store the order of PDBs in the dendrogram
dendrogram_order <- data.frame(pdb_id = pdb_id_labels,
                      order = seq(1,length(pdb_id_labels),1))

# Saving Dendrogram Order
write.csv(dendrogram_order, file = file.path(data_path, "dendrogram_order.csv"), row.names = FALSE)


if (is.na(custom_cut_height) || custom_cut_height <= 0 || is.null(custom_cut_height)) {
  
  # Setting cut_height to the mean distance
  cut_height <- mean(cluster_dist)
  print(paste0("Using mean Euclidean distance as the cut height: ", cut_height))
  
} else {
  
  cut_height = custom_cut_height
  print(paste0("Using custom cut height: ", cut_height))

}

# Extracting pdb_id from each cluster group
clusters <- cutree(hc_average, h = cut_height)
grouped_pdbs <- split(df_wide$pdb_id, clusters)

# Convert the list to a dataframe
cluster_groups <- do.call(rbind, lapply(names(grouped_pdbs), function(name) {
  data.frame(group = name, pdb_id = grouped_pdbs[[name]], stringsAsFactors = FALSE)
}))

# Ordering the colour column by the RMSD Dendrogram order
pdb_id_labels_factor <- factor(pdb_id_labels, levels = pdb_id_labels)
cluster_groups <- cluster_groups[order(factor(cluster_groups$pdb_id, levels = levels(pdb_id_labels_factor))), ]

# Ordering cluster_groups by dendrogram order 
cluster_groups <- cluster_groups[match(pdb_id_labels, cluster_groups$pdb_id), ]

tmp <- as.character(unique(cluster_groups$group))

cluster_groups_sorted <- cluster_groups

for (i in 1:length(tmp)) {
  
  cluster_groups_sorted$group[cluster_groups_sorted$group == tmp[i]] <- paste0(i,".")  
  
}

# Removing temporary "."
cluster_groups_sorted$group <- gsub("\\.", "", cluster_groups_sorted$group)


# Saving Cluster Groups
write.csv(cluster_groups_sorted, file = file.path(data_path, "RMSD_cluster_groups.csv"), row.names = FALSE)

# Making group column numeric
cluster_groups$group <- as.numeric(cluster_groups$group)

# Setting branch colours 
colours <- c("red", "blue", "gold2", "purple", "darkorange2", "cyan", "green", "pink", "chocolate4", 
             "darkgreen", "magenta", "deepskyblue", "darkviolet", "coral", "darkslategray")

# Extend the color vector to match groups
branch_colours <- rep(colours, length.out = as.numeric(max(cluster_groups$group)))

# Using group number to index the colour palette
cluster_groups$colour <- branch_colours[as.numeric(cluster_groups$group)]

# Taking only the unique colour values as the plot requires one value for each group (not for each line)
group_colours <- unique(cluster_groups$colour)


# Plotting
avg_dend_obj <- as.dendrogram(hc_average)
avg_col_dend <- color_branches(avg_dend_obj, h = cut_height, col = group_colours)
labels(avg_col_dend) <- pdb_id_labels
avg_col_dend <- set(avg_col_dend, "branches_lwd", 2.5) %>% set("labels_cex", 1)

# Save the cluster plot as a PNG file
png(file.path(output_dir, "RMSD_cluster_dendrogram.png"), width = 6000, height = 1500, res = 300)
plot(avg_col_dend, ylab = "Euclidean Distance", xlab = "PDB ID",
     cex.lab = 1.2, font.lab = 2, cex.axis = 1)
abline(h = cut_height, col = "grey30", lty = 2, lw = 2)
title("Cluster Dendrogram - RMSD")
dev.off()

#############################
### RMSD Heatmap Plotting ###
#############################

# RMSD = Root Mean Squared Distance

# Converting dataframe into long format for plotting
mean_distance_heatmap <- df_wide %>% pivot_longer(!pdb_id, names_to = "comparison", values_to = "mean_distance")

# Ordering by cluster order
mean_distance_heatmap$pdb_id <- factor(mean_distance_heatmap$pdb_id, levels = pdb_id_labels)
mean_distance_heatmap$comparison <- factor(mean_distance_heatmap$comparison, levels = pdb_id_labels)

# Heatmap plot
ggplot(mean_distance_heatmap, aes(x = comparison, y = pdb_id, fill= mean_distance)) +
  geom_tile() +
  ggtitle("RMSD Heatmap") +
  ylab("PDB") +
  xlab("PDB") +
  #scale_fill_gradient(low="black", high="red", na.value = "white") +
  scale_fill_gradientn(colors = c("red", "black", "blue")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle = 90,
                                   vjust = 0.4, hjust = 1),
        axis.title=element_text(size = 25, face = "bold"),
        axis.line.x = element_line(color="black", linewidth = 1.5),
        axis.line.y = element_line(color="black", linewidth = 1.5),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 16, face = "bold")) +
  labs(fill = "RMSD")

# Saving the heatmap
ggsave(plot = last_plot(), file = file.path(output_dir, "RMSD_heatmap.png"),
      width = 18, height = 15)

# Saving data
write.csv(mean_distance_heatmap, file = file.path(data_path, "RMSD_heatmap.csv"), row.names = FALSE)


# Cleaning Environment
#rm(list = ls())
invisible(gc())


