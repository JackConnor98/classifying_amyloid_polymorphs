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

# Convert 'plot' parameter to logical
plot <- tolower(config$plot)  # Ensure case insensitivity

# Convert 'custom_cut_height' parameter to numeric
custom_cut_height <- as.numeric(config$custom_cut_height)

# Check if the configuration was correctly read
if (plot == "TRUE") {
  plot <- TRUE
  print("Plotting is enabled.")
} else {
  plot <- FALSE
  print("Plotting is disabled.")
}

# Check the custom_cut_height
print(paste("Custom cut height:", custom_cut_height))

############################################################################################

# Creating Directories
if(!dir.exists("Output/RMSD")) dir.create("Output/RMSD")
if(!dir.exists("Output/RMSD/single_reference")) dir.create("Output/RMSD/single_reference")

############################################################################################

# Read in df
df <- read.csv("Output/PDBs/Alignments/all_unique_chains_aligned.csv")

# Remove duplicated rows based on the specified columns (occurs once at residue 42 in 7wmm)
df <- df[!duplicated(df[, c('pdb_id', 'pdb', 'chain', 'atom_id', 'residue_number')]), ] 
# removed 'residue_name' as this breaks analysis for some Tau structures e.g. 9eo7 which is a V337M mutant but both V and M are present at 337 in the PDB...

############################################################################################

### Removing regions outside the fibril core ###

metadata <- read_delim("Output/selected_pdbs_metadata.txt",
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)

metadata <- subset(metadata, select = c("PDB ID",  "Residues Ordered"))
names(metadata) <- c("pdb", "ordered_residues")

metadata$core_start <- NA
metadata$core_end <- NA

for (i in 1:nrow(metadata)) {

  range <- metadata$ordered_residues[i]
  range <- str_split_1(range, paste("-", ",", sep = "|"))

  metadata$core_start[i] <- as.numeric(range[1])
  metadata$core_end[i] <- as.numeric(range[length(range)])

}

# Combining df and metadata
df <- merge(df, metadata, by = "pdb")

# Changing Start and End NAs for non-amyloid atlas pdbs
df <- df %>%
  mutate(core_start = replace_na(core_start, 0),
         core_end = replace_na(core_end, Inf))

# Filtering residues that are outside the core
df <- df %>% filter(residue_number >= core_start,
                    residue_number <= core_end)

############################################################################################

### Removing Low Resolution Residues ###

# Importing data
high_resolution_residues <- read.csv("Output/Validation/high_resolution_residues.csv")

# Renaming columns to match df
names(high_resolution_residues)[names(high_resolution_residues) == "fibril"] <- "fibril_number"
names(high_resolution_residues)[names(high_resolution_residues) == "resno"] <- "residue_number"
names(high_resolution_residues)[names(high_resolution_residues) == "pdb_id"] <- "pdb_id"

# Merging with main data frame
df <- merge(df, high_resolution_residues, by = c("pdb_id", "fibril_number", "residue_number"))

############################################################################################

# Inter-PDB RMSD

cat("Calculating Inter-PDB RMSD\n")

# Getting amyloid names from first chain only
amyloid_names <- unique(df$pdb_id)

#Initializing data frame to store mean_distance scores for each comparison
mean_distance_data <- data.frame(pdb_id = unique(df$pdb_id))


for (i in amyloid_names) {

  # Remove any positions where reference chain is not found (Unable to compare residues that don't exist in reference chain)
  filtered_df <- df %>%
                  group_by(residue_number) %>%
                    filter(any(pdb_id == i))

  # Creates new columns that equals the coordinates of the reference chain for each pdb
  filtered_df <- filtered_df %>%
    group_by(residue_number) %>%
    mutate(x_comp = x[pdb_id == i],
           y_comp = y[pdb_id == i],
           z_comp = z[pdb_id == i])

  # Finding distance between two points
  filtered_df <- filtered_df %>%
                  mutate(distance = sqrt((x_comp - x)^2 + (y_comp - y)^2 + (z_comp - z)^2))

  if (plot) {

    # Plotting distance

    ggplot(filtered_df, aes(x = residue_number, y = distance, colour = pdb_id)) +
      geom_point(size = 3, alpha = 0.5) +
      geom_line(aes(group = interaction(pdb_id, cumsum(c(0, diff(residue_number) != 1)))),
                linewidth = 1, alpha = 0.5) +
      ggtitle(paste(i, " - Chain Comparison", sep = "")) +
      ylab("Distance") +
      xlab("Residue Position") +
      scale_x_continuous(breaks = seq(0,140,2)) +
      theme(panel.grid.major = element_line(colour = "grey30", linewidth = 0.5,
                                            linetype = "dashed"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
            plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
            axis.text.y = element_text(size = 16, colour = "black"),
            axis.text.x = element_text(size = 16, colour = "black"),
            axis.title=element_text(size = 25, face = "bold"),
            axis.line.x = element_line(color="black", linewidth = 1.5),
            axis.line.y = element_line(color="black", linewidth = 1.5),
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 10, face = "bold"))

    ggsave(plot = last_plot(), file = paste("Output/RMSD/single_reference/", i, "_distance.png", sep = ""),
           width = 20, height = 8)
    
  }

  #################################################################################################
  #################################################################################################
  #################################################################################################
  #################################################################################################

  # Creating data frame for mean_distance heatmap

  # Finding the mean distance (i.e mean_distance) for each pdb
  comp_mean_distance <- filtered_df %>% group_by(pdb_id) %>% summarise_at(vars(distance), funs(mean))


  if (plot) {
    # Plotting comp_mean_distance

    ggplot(comp_mean_distance, aes(x = reorder(pdb_id, distance), y = distance)) +
      geom_point(size = 4) +
      ggtitle(paste(i, " - Chain Comparison", sep = "")) +
      ylab("RMSD") +
      xlab("PDB") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(comp_mean_distance$distance) * 1.05),
                         breaks = seq(0,1000,5)) +
      theme(panel.grid.major = element_line(colour = "grey30", linewidth = 0.5,
                                            linetype = "dashed"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
            plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
            axis.title = element_text(size = 18, face = "bold", colour = "black"),
            axis.text.x = element_text(size = 14, colour = "black",
                                       angle = 90, vjust = 0.4, hjust = 1),
            axis.text.y = element_text(size = 14, colour = "black"))

    ggsave(plot = last_plot(), file = paste0("Output/RMSD/single_reference/", i, "_RMSD.png"), width = 14, height = 8)

  }

  # Adding mean_distance to data frame
  mean_distance_data <- merge(mean_distance_data, comp_mean_distance, all = TRUE)

  # renaming last column so I know what comparison the similarities refer to
  names(mean_distance_data)[length(names(mean_distance_data))] <- i

}

# Saving Mean Chain Distance Data
write.csv(mean_distance_data, "Output/RMSD/all_RMSD_comparisons.csv", row.names=FALSE)

################################
################################
################################

### Clustering Based On RMSD ###

################################
################################
################################

# Removing reference columns
cluster_data <- subset(mean_distance_data, select = setdiff(names(mean_distance_data), c("pdb_id")))

# Setting NA values to the max distance - This occurs when two PDBs have no overlapping residues 
#NA_distance <- max(cluster_data, na.rm = TRUE)

# Setting NA values to the mean + 3SD - This occurs when two PDBs have no overlapping residues 
NA_distance <- mean(as.matrix(cluster_data), na.rm = TRUE) + 3*sd(as.matrix(cluster_data), na.rm = TRUE)

cluster_data[is.na(cluster_data)] <- NA_distance

# Scaling Cluster Data (makes mean of all columns = 0 and SD = 0-1)
scaled_cluster_data <- as.data.frame(scale(cluster_data))
#summary(scaled_cluster_data)

# Calculating Distance
cluster_dist <- dist(scaled_cluster_data, method = "euclidean")

# Clustering
hc_average <- hclust(cluster_dist, method = "average")


### Scree Plot to assigning cut height ###
scree_data <- hc_average$height %>%
              as_tibble() %>%
              add_column(groups = length(hc_average$height):1) %>%
              rename(height=value)

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
ggsave(plot = last_plot(), filename = "Output/RMSD/Scree_plot.png", height = 4, width = 8)



# Making labels
pdb_id_labels <- mean_distance_data$pdb_id[hc_average$order]

# Creating a dataframe to store the order of PDBs in the dendrogram
dendrogram_order <- data.frame(pdb_id = pdb_id_labels,
                      order = seq(1,length(pdb_id_labels),1))

# Saving Dendrogram Order
write.csv(dendrogram_order, file = "Output/RMSD/dendrogram_order.csv", row.names = FALSE)


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
grouped_pdbs <- split(mean_distance_data$pdb_id, clusters)

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
write.csv(cluster_groups_sorted, file = "Output/RMSD/RMSD_cluster_groups.csv", row.names = FALSE)

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
png("Output/RMSD/RMSD_cluster_dendrogram.png", width = 5000, height = 1500, res = 300)
plot(avg_col_dend, ylab = "Euclidean Distance", xlab = "PDB ID",
     cex.lab = 1.2, font.lab = 2, cex.axis = 1)
abline(h = cut_height, col = "grey30", lty = 2, lw = 2)
#title("Cluster Dendrogram - RMSD")
dev.off()

#############################
### RMSD Heatmap Plotting ###
#############################

# RMSD = Root Mean Squared Distance

# Converting dataframe into long format for plotting
mean_distance_heatmap <- mean_distance_data %>% pivot_longer(!pdb_id, names_to = "comparison", values_to = "mean_distance")

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
ggsave(plot = last_plot(), file = paste("Output/RMSD/RMSD_heatmap.png", sep = ""),
       width = 18, height = 15)

# Saving data
write.csv(mean_distance_heatmap, file = "Output/RMSD/RMSD_heatmap.csv", row.names = FALSE)


# Cleaning Environment
#rm(list = ls())
invisible(gc())


