cat("Running: defining_stable_regions.R\n")

# Installing and Loading Packages
required_packages <- c("dplyr", "stringr", "readr", "ggplot2", "tidyr", "zoo", "pastecs")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]   

# Install missing packages
if(length(not_installed)) install.packages(not_installed)    

# Load all required packages
invisible(lapply(required_packages, require, character.only = TRUE))

# Cleaning Environment
rm(list = ls())

################################################################################

# Creating Directories
if(!dir.exists(file.path("Output", "stable_regions"))) dir.create(file.path("Output", "stable_regions"))

################################################################################

# Importing Data
filtered_df <- read.csv(file.path("Output", "thermodynamics", "foldx_stability_filtered.csv"))

# Importing the mean delta G per residue for each PDB 
PDB_mean_per_res <- read.csv(file.path("Output", "thermodynamics", "mean_deltaG_per_residue.csv"))

# Playing with different sliding windows

window_size <- 3

PDB_mean_per_res <- PDB_mean_per_res %>%
                      group_by(pdb_id) %>%
                        mutate(sliding_window = rollapply(mean_energy, width = window_size, FUN = mean, align = "center", fill = NA)) %>% 
                          ungroup()

################################################################################

# Finding local minima (troughs)
mean_at_each_pos <- PDB_mean_per_res %>% group_by(Pos) %>% summarise(mean_Pos_sliding_window = mean(sliding_window, na.rm = TRUE))
mean_at_each_pos <- na.omit(mean_at_each_pos)


mean_stability <- mean(PDB_mean_per_res$sliding_window, na.rm = TRUE)

# Get the turning points object
tp <- pastecs::turnpoints(mean_at_each_pos$mean_Pos_sliding_window, calc.proba = TRUE)

# Creating a turnpoint df from tp
info <- extract(tp, peak = 1, pit = -1)
type <- info[info != 0]

tp_df <- data.frame(Pos = mean_at_each_pos$Pos[tp$tppos], # tppos counts each data point, since my data doesnt start at 1 I need to account for it
                    info = as.character(type),
                    probability = tp$proba)

# Define a probability threshold
threshold <- 0.15  # for instance, if you're interested in low probability values

# Filter the turning points object based on probability
tp_df <- tp_df %>% filter(probability < threshold)

# Create a data frame for plotting
mean_at_each_pos <- merge(mean_at_each_pos, tp_df, by = "Pos", all = TRUE)
mean_at_each_pos$info[is.na(mean_at_each_pos$info)] <- "0"

#####################################################################################################

### Selecting Stable Regions ###

# Selecting peaks above the mean line
significant_peaks <- mean_at_each_pos$Pos[mean_at_each_pos$info == 1]

# Initialize variables
stabilisers <- c()
region_list <- c()
region <- 1
stabiliser_count <- 0

# Loop through the data to find stabilizers
for (i in seq_along(mean_at_each_pos$mean_Pos_sliding_window)) {
  
  # Check if the current position is stabilizing
  if (mean_at_each_pos$mean_Pos_sliding_window[i] < mean_stability & mean_at_each_pos$info[i] != 1) {
    stabilisers <- c(stabilisers, mean_at_each_pos$Pos[i])
    region_list <- c(region_list, region)
    stabiliser_count <- stabiliser_count + 1
  }
  
  # If a significant peak is hit, check if the current region has more than 1 stabilizer
  if (mean_at_each_pos$Pos[i] %in% significant_peaks) {
    
    # Only proceed if there are more than 1 stabilizers in the region
    if (stabiliser_count > 1) {
      region <- region + 1
    } else {
      # Remove stabilisers with less than 2 occurrences
      stabilisers <- stabilisers[region_list != region]
      region_list <- region_list[region_list != region]
    }
    
    # Reset stabilizer count for the next region
    stabiliser_count <- 0
  }
}

# Creating a data frame of only significantly stabilising positions and which group
stabilising_positions <- data.frame(Pos = stabilisers,
                          region = region_list)

# Fixing the group numbers to start from 1
unique_regions <- unique(region_list)

for (i in seq_along(unique_regions)) {
  
  stabilising_positions$region[stabilising_positions$region == unique_regions[i]] <- i
  
}

# Summarising into stable regions
stable_regions <- stabilising_positions %>% 
                    group_by(region) %>% 
                      summarise(min = min(Pos), max = max(Pos))

stable_regions <- stable_regions %>% mutate(residues = paste0(min, "-", max))

stable_regions$region <- as.character(stable_regions$region)

stable_regions$region <- factor(stable_regions$region, levels = sort(as.numeric(unique(stable_regions$region))))

# Assigning colours

colours <- c("red", "blue", "gold2", "purple", "darkorange2", "cyan", "green", "pink", "chocolate4", 
                    "darkgreen", "magenta", "deepskyblue", "darkviolet", "coral", "darkslategray")

# Extend the color vector to match `region`
region_colours <- rep(colours, length.out = region)


# Visualising Stabilising Regions
ggplot(mean_at_each_pos, aes(x = Pos, y = mean_Pos_sliding_window)) +
  geom_rect(data = stable_regions, aes(xmin = min - 0.5, xmax = max + 0.5, ymin = -Inf, ymax = Inf, fill = region),
            alpha = 0.75, inherit.aes = FALSE)+#, fill = "grey50") +
  geom_line(color = "black", linewidth = 1) +
  geom_point(aes(color = info), size = 2) +
  scale_fill_manual(values = region_colours) +
  scale_color_manual(values = c("1" = "blue", "0" = "black", "-1" = "red"),
                     label= c("Trough", "NA", "Peak")) +
  scale_x_continuous(breaks = seq(0,1000,10)) +
  geom_hline(yintercept = mean_stability, linetype = "solid", size = 1) +  # Mean line
  labs(title = "Turning Points", x = "Residue Number", y = "Mean Pos Sliding Window", colour = "Type", fill = "Region") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        plot.title = element_text(size = 14, colour = "black", face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 10, colour = "black")) 

ggsave(plot = last_plot(), file = file.path("Output", "stable_regions", "turning_points_and_regions.png"), height = 5, width = 8, bg = "white")

ggplot(data = PDB_mean_per_res, aes(x=Pos)) +
  geom_rect(data = stable_regions, aes(xmin = min - 0.5, xmax = max + 0.5, ymin = -Inf, ymax = Inf, fill = region),
            alpha = 0.75, inherit.aes = FALSE) +
  scale_fill_manual(values = region_colours) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 1, linetype = "solid", alpha = 0.9) +
  geom_line(aes(y = sliding_window, group = interaction(pdb_id, cumsum(c(0, diff(Pos) != 1)))), 
                 linewidth = 0.5, alpha = 0.25) +
  ggtitle(paste0("Sliding Window of ", window_size)) +
  ylab(bquote(bold(Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
  xlab("Residue Number") +
  scale_x_continuous(breaks = seq(0,99999, by = 10)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        plot.title = element_text(size = 14, colour = "black", face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none")

ggsave(plot = last_plot(), file = file.path("Output", "stable_regions", "sliding_window_and_regions.png"), height = 4, width = 8, bg = "white")


ggplot(data = PDB_mean_per_res, aes(x=Pos)) +
  geom_rect(data = stable_regions, aes(xmin = min - 0.5, xmax = max + 0.5, ymin = -Inf, ymax = Inf, fill = region),
            alpha = 0.75, inherit.aes = FALSE) +
  scale_fill_manual(values = region_colours) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 1, linetype = "solid", alpha = 0.9) +
  geom_line(aes(y = mean_energy, group = interaction(pdb_id, cumsum(c(0, diff(Pos) != 1)))), 
            linewidth = 0.5, alpha = 0.25) +
  ylab(bquote(bold(Delta * "G"*degree * " per residue (kcal.mol"^-1 * ")"))) +
  xlab("Residue Number") +
  scale_x_continuous(breaks = seq(0,99999, by = 10)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none")

ggsave(plot = last_plot(), file = file.path("Output", "stable_regions", "raw_delta_G_and_regions.png"), height = 4, width = 8, bg = "white")

# Saving stable regions
stable_regions_minimal <- stable_regions %>% select("residues", "region")
write.csv(stable_regions_minimal, file = file.path("Output", "stable_regions", "stable_regions.csv"), row.names = FALSE)
