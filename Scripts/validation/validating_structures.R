cat("Running: validating_structures.R\n")

# Script Logic:
  # Q-scores for each PDB are read and combined into q_score_df

  # The mean and SD for Q-score is calculated using all residues from every PDB
  # The Q-score threshold is set to mean - 1SD
  # The average Q-score for each PDB is calculated
  # PDBs with a mean Q-score below threshold are removed
  # Mean Q-score is recalculated from the remaining PDBs
  # For the remaining PDBs: residues with a Q-score below threshold are removed

  # Issue - Chains within the same fibril often have some residues found in some and not others due to filtering
  # Current solution - if a residue is of good resolution for at least 1 chain in the fibril - keep it

  # Create a dataframe containing the pdb_id, fibril number and residue number for those that are of good resolution

############################################################################################

# Installing and Loading Packages
required_packages <- c("dplyr", "stringr", "tidyr", "readr", "ggplot2", "XML", "xml2", "bio3d")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]   

# Install missing packages
if(length(not_installed)) install.packages(not_installed)    

# Load all required packages
invisible(lapply(required_packages, require, character.only = TRUE))

# Cleaning Environment
rm(list = ls())

############################################################################################

# Locate all PDB Q-score data in folder that need to be looped through
folder_path <- file.path("Output", "Validation")

data_path <- paste0(folder_path, "/data")

# Get a list of all .csv files in the folder
xml_files <- list.files(path = data_path, pattern = "\\.gz$", full.names = TRUE)

############################################################################################

save_path <- file.path(folder_path, "plots")
single_pdb_path <- file.path(save_path, "single_pdbs")

# Creating Directories
if(!dir.exists(save_path)) dir.create(save_path)
if(!dir.exists(single_pdb_path)) dir.create(single_pdb_path)

############################################################################################

# Importing metadata
selected_pdbs_metadata <- read_delim(file.path("Output", "selected_pdbs_metadata.txt"), 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
# Fixing names
names(selected_pdbs_metadata)[names(selected_pdbs_metadata) == "PDB ID"] <- "PDB"
names(selected_pdbs_metadata)[names(selected_pdbs_metadata) == "Resol- ution (Å)"] <- "resolution"

# Selecting columns
resolution_df <- selected_pdbs_metadata %>% select(c("PDB", "resolution"))

############################################################################################

##################################################################################
### Creating a dataframe containing the Q-scores for every residue of each PDB ###
##################################################################################

validation_data <- data.frame()

for (file in xml_files) {
  
  # Getting filename
  filename <- tail(str_split_1(file, "/"), 1)
  
  # Read the XML file
  tmp <- read_xml(file)
  
  # Find all ModelledSubgroup nodes
  nodes <- xml_find_all(tmp, ".//ModelledSubgroup")
  
  # Get all unique attribute names
  all_attributes <- unique(unlist(lapply(nodes, function(node) names(xml_attrs(node)))))
  
  # Checking to see if Validation data contains Q-scores
  if ("Q_score" %in% all_attributes) {
    
    # Extract attributes and ensure consistency
    df <- do.call(rbind, lapply(nodes, function(node) {
      attrs <- xml_attrs(node)
      # Ensure all columns are present
      df_row <- as.data.frame(t(attrs), stringsAsFactors = FALSE)
      missing_cols <- setdiff(all_attributes, names(df_row))
      df_row[missing_cols] <- NA
      df_row
    }))
    
    # Convert to a data frame
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    
    # Adding pdb column
    df$pdb <- str_split_1(filename, "_")[1]
    
    
    # Merging to main df
    validation_data <- bind_rows(validation_data, df)
    
  } else {
    
    print(paste0(str_split_1(filename, "_")[1], " : No Q-scores found"))
    
  }
  
  
}

# Selecting key columns 
q_score_df <- validation_data %>% subset(select = c("pdb", "Q_score", "chain", "resnum", "resname"))

# Changing resnum to resno to better fit other data
names(q_score_df)[names(q_score_df) == "pdb"] <- "PDB"
names(q_score_df)[names(q_score_df) == "resnum"] <- "resno"

# Converting columns to numeric
q_score_df[, c("Q_score", "resno")] <- lapply(q_score_df[, c("Q_score", "resno")], as.numeric)

# Removing NAs (This will remove ssNMR structures as Q-scores are for cryoEM)
q_score_df <- q_score_df[!is.na(q_score_df$Q_score), ]

# Character vector containing the capitalized 3-letter codes for all amino acids in all capitals
amino_acids <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", 
                 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")

# Removing non-amino acid components
q_score_df <- q_score_df %>% filter(resname %in% amino_acids)

### Adding in pdb_id ###

# Importing data
pdb_info <- read.csv(file.path("Output", "PDBs", "COM_and_fibril.csv"))

# Selecting columns of interes
pdb_info <- pdb_info %>% select(c("pdb_id", "PDB", "chain", "fibril"))

# Adding to q_score_df
q_score_df <- merge(pdb_info, q_score_df, by = c("PDB", "chain"))

# Saving Q-Score data frame
write.csv(q_score_df, file.path(folder_path, "Q_Scores.csv"), row.names = FALSE)

##############################################
### Filtering Out Poorly Resolved Residues ###
##############################################

# Merging resolution with q_scores
q_score_df <- merge(q_score_df, resolution_df, by = "PDB")


### Optional ###
# Setting a hard limit for reported resolution to 4.0A #
#q_score_df <- q_score_df %>% filter(resolution < 4)


# Calculating the mean and standard deviation for the Q-score across all structures
# Threshold for Q-score will me set to mean - 1SD

# Calculating mean and standard deviation
mean_Q_score <- mean(q_score_df$Q_score, na.rm = TRUE)
sd_Q_score <- sd(q_score_df$Q_score, na.rm = TRUE)
Q_score_threshold <- mean_Q_score - sd_Q_score


####################################################
### Calculating the mean Q-score for each pdb_id ###
####################################################

PDB_q_score <- q_score_df %>%
                group_by(pdb_id) %>%
                  summarise(mean_Q_score = mean(Q_score))

# Putting in descending order
PDB_q_score <- PDB_q_score[order(PDB_q_score$mean_Q_score, decreasing = TRUE),]

# Saving the mean Q-score for each PDB
write.csv(PDB_q_score, file = file.path(folder_path, "mean_PDB_Q_score.csv"), row.names = FALSE)

### Plotting mean PDB Q-score ###

# Ordering 'pdb' by 'mean_Q_score'
PDB_q_score$pdb_id <- factor(PDB_q_score$pdb_id, levels = PDB_q_score$pdb_id[order(PDB_q_score$mean_Q_score, decreasing = TRUE)])

# Plotting
ggplot(PDB_q_score, aes(x = pdb_id, y = mean_Q_score)) +
  geom_hline(yintercept = mean_Q_score, colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = Q_score_threshold, colour = "red", linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 1) +
  labs(title = "",
       x = "PDB",
       y = "Mean Q-Score") +
  theme_classic() +
  theme(panel.border = element_rect(linewidth = 1, fill = NA),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey50", linetype = "dashed"),
        axis.title = element_text(size = 18, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black", 
                                   angle = 90, vjust = 0.25))

ggsave(file = file.path(save_path, "mean_Q_score.png"), plot = last_plot(), height = 4, width = 8, bg = "white")

###############################################################
### Plotting each PDB individually for quality/sense checks ###
###############################################################

pdb_names <- unique(q_score_df$pdb_id)

max_y <- max(q_score_df$Q_score)
min_y <- min(q_score_df$Q_score)
max_x <- max(q_score_df$resno)
min_x <- min(q_score_df$resno)

for (i in pdb_names) {
  
  # Filtering current PDB
  x <- q_score_df %>% filter(pdb_id == i)
  
  # Colouring outliers
  x$q_outliers <- ifelse(x$Q_score < Q_score_threshold, "outlier", "good")
  
  # Getting Mean Q-score
  current_mean_q_score <- PDB_q_score$mean_Q_score[PDB_q_score$pdb_id == i]
  current_mean_q_score <- signif(current_mean_q_score, digits = 3)
  
  # Getting Resolution
  current_resolution <- unique(q_score_df$resolution[q_score_df$PDB == i])
  
  # Plotting Q-Score
  ggplot() +
    geom_point(data = x, aes(x = resno, y = Q_score, colour = q_outliers), size = 1) +
    geom_hline(yintercept = mean_Q_score, linewidth = 1, linetype = "solid", color = "black") +
    geom_hline(yintercept = mean_Q_score - sd_Q_score, linewidth = 1, linetype = "dashed", colour = "grey30") +
    scale_colour_manual(values = c("outlier" = "grey70", "good" = "black")) +
    labs(title = paste0("PDB = ", i, 
                        "  |  Resolution = ", current_resolution, "\u00C5", 
                        "  |  Mean Q-Score = ", current_mean_q_score),
         x = "Residue Number",
         y = "Q Score") +
    xlim(c(min_x, max_x)) +
    ylim(c(min_y, max_y)) +
    theme_classic() +
    theme(panel.border = element_rect(linewidth = 1, fill = NA),
          plot.title = element_text(size = 16, hjust = 0.5),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 14, colour = "black"),
          legend.position = "none")
  
  ggsave(file = file.path(single_pdb_path, paste0(i, "_Q_score.png")), plot = last_plot(), height = 4, width = 8, bg = "white")
  
}

######################################
### Removing Low Mean Q-score PDBs ###
######################################

# Selecting PDBs with a mean Q-score > Q_score_threshold
PDB_q_score$pdb_id <- as.character(PDB_q_score$pdb_id)
high_resolution_PDBs <- PDB_q_score$pdb_id[PDB_q_score$mean_Q_score > Q_score_threshold]

# Saving a list of high_resolution PDBs
high_resolution_PDBs_df <- PDB_q_score[PDB_q_score$mean_Q_score > Q_score_threshold,]
write.table(high_resolution_PDBs_df, file = file.path("Output", "Validation", "high_resolution_pdb_ids.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)


filtered_df <- q_score_df %>% filter(pdb_id %in% high_resolution_PDBs)

filtered_df <- filtered_df %>% select("pdb_id", "chain", "resno", "Q_score")


# Recalculating mean and standard deviation with low resolution PDBs removed
mean_Q_score <- mean(filtered_df$Q_score, na.rm = TRUE)
sd_Q_score <- sd(filtered_df$Q_score, na.rm = TRUE)
Q_score_threshold <- mean_Q_score - sd_Q_score

#####################################
### Removing Low Q-score Residues ###
#####################################

# Selecting good resolution residues
good_resolution <- filtered_df %>% filter(Q_score >= mean_Q_score - sd_Q_score)

################################################
### Visualising Q-score Vairance in Residues ###
################################################

### Density Plot ###

# Calculate the density
density_data <- density(filtered_df$Q_score)

# Create a dataframe from density data
density_df <- data.frame(x = density_data$x, y = density_data$y)

# Saving density data
write.csv(density_df, file.path("Output", "Validation", "q_score_density_data.csv"), row.names = FALSE)


# Plot with different colors on either side of the dashed line
ggplot() +
  geom_area(data = subset(density_df, x <= Q_score_threshold), 
            aes(x = x, y = y), 
            fill = "grey50", alpha = 0.75, 
            colour = "black", linewidth = 0.8) +
  geom_area(data = subset(density_df, x > Q_score_threshold), 
            aes(x = x, y = y), 
            fill = "blue", alpha = 0.75, 
            colour = "black", linewidth = 0.8) +
  geom_vline(xintercept = mean_Q_score, colour = "black", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = Q_score_threshold, colour = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Q Score",
       y = "Density") +
  scale_x_continuous(breaks = seq(-1,1,0.2), expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(panel.border = element_rect(linewidth =  1, fill = NA),
        axis.title = element_text(size = 18, face = "bold", colour = "black"),
        axis.text = element_text(size = 14, colour = "black"))

ggsave(file = file.path((save_path, "q_score_density.png"), plot = last_plot(), height = 4, width = 6, bg = "white")

### Stacked Bar Chart ###

# Counting the number of occurrences for each residue
good_resolution_counts <- good_resolution %>% count(resno, name = "good_count")
total_counts <- filtered_df %>% count(resno, name = "total_count")

# Merging good and total counts
count_df <- merge(good_resolution_counts, total_counts, by = "resno")

# Calculating % good
count_df <- count_df %>% mutate(percentage = good_count / total_count * 100,
                                bad_count = total_count - good_count)

# Transform to long format
long_count_df <- count_df %>%
  pivot_longer(cols = c("good_count", "bad_count"), 
               names_to = "count_type", 
               values_to = "count")

# Saving bar chart data
write.csv(long_count_df, file.path("Output", "Validation", "resolution_bar_chart_data.csv"), row.names = FALSE)

# Create the stacked bar plot
ggplot(long_count_df, aes(x = factor(resno), y = count, fill = count_type)) +
  geom_bar(stat = "identity", alpha = 0.75, colour = "black") +
  scale_fill_manual(values = c("bad_count" = "grey", "good_count" = "blue"),
                    labels = c("Discarded", "Accepted")) +
  labs(title = "",
       x = "Residue Number",
       y = "Count",
       fill = "Resolution") +
  scale_x_discrete(breaks = seq(0,max(long_count_df$resno), 10)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(panel.border = element_rect(linewidth =  1, fill = NA),
        axis.title = element_text(size = 18, face = "bold", colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, face = "bold", colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),)

ggsave(file = file.path(save_path, "residue_count_stacked_bar.png"), plot = last_plot(), height = 4, width = 8, bg = "white")


# Handling chains with different number of residues
# Current solution - if the residue occurs in at good resolution at least once in a fibril - it is included

# Adding back in other information
pdb_info_no_chain <- pdb_info %>% select(-c("chain")) %>% distinct()
tmp <- merge(pdb_info_no_chain, good_resolution, by = "pdb_id")

# Selecting resno's that appear at least once in each fibril
high_resolution_residues <- tmp %>% 
  group_by(pdb_id, fibril) %>%
  distinct(resno)

# Putting in ascending resno order
high_resolution_residues <- high_resolution_residues[order(high_resolution_residues$pdb_id, high_resolution_residues$fibril, high_resolution_residues$resno), ]

# Saving final data
write.csv(high_resolution_residues, file = file.path(folder_path, "high_resolution_residues.csv"), row.names = FALSE)











