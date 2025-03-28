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
if(!dir.exists("Output/stable_regions")) dir.create("Output/stable_regions")
if(!dir.exists("Output/stable_regions/normality")) dir.create("Output/stable_regions/normality")

################################################################################

# Importing Data
filtered_df <- read.csv("Output/thermodynamics/foldx_stability_filtered.csv")

# Importing the mean delta G per residue for each PDB 
PDB_mean_per_res <- read.csv("Output/thermodynamics/foldx_stability_filtered_normalised.csv")

# Playing with different sliding windows

window_size <- 5

PDB_mean_per_res <- PDB_mean_per_res %>%
                      group_by(pdb_id) %>%
                        mutate(sliding_window = rollapply(mean_energy, width = window_size, FUN = mean, align = "center", fill = NA)) %>% 
                          ungroup()

################################################################################

### Normality Testing ###

# Histogram
ggplot(PDB_mean_per_res, aes(x = sliding_window)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
  labs(title = paste0("Sliding Window of ", window_size, " Averaged"), x = "Mean \u0394G per residue (kcal.mol-1)", y = "Density") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none")

# Saving plot
ggsave(plot = last_plot(), file = "Output/stable_regions/normality/histogram.png", height = 4, width = 6)


# Q-Q plot

# Open a PNG graphics device
png(filename = "Output/stable_regions/normality/QQ_plot.png",
    width = 6, height = 4, units = "in", res = 600)  

qqnorm(PDB_mean_per_res$sliding_window)
qqline(PDB_mean_per_res$sliding_window, col = "red")

# Close the PNG graphics device
dev.off()

# Saving normality test results
sink("Output/stable_regions/normality/kolmogorov_smirnov_test.txt")

# Normality test (Kolmogorov-Smirnov)
ks.test(PDB_mean_per_res$sliding_window, "pnorm",
        mean = mean(PDB_mean_per_res$sliding_window, na.rm = TRUE),
        sd = sd(PDB_mean_per_res$sliding_window, na.rm = TRUE))

sink()

################################################################################

# Function to find ranges with gaps <= 2  
# find_stable_regions <- function(x) {
#   if (length(x) == 0) return(character(0))
#   x <- sort(x)
#   ranges <- c()
#   start <- x[1]
#   end <- x[1]
#   
#   for (i in 2:length(x)) {
#     if (x[i] <= end + 2) {
#       end <- x[i]
#     } else {
#       if (start == end) {
#         ranges <- c(ranges, as.character(start))
#       } else {
#         ranges <- c(ranges, paste0(start, "-", end))
#       }
#       start <- x[i]
#       end <- x[i]
#     }
#   }
#   
#   if (start == end) {
#     ranges <- c(ranges, as.character(start))
#   } else {
#     ranges <- c(ranges, paste0(start, "-", end))
#   }
#   
#   return(ranges)
# }

################################################################################


# Using mean stability for each residue for all pdb_ids (to account for different number of chains)

# Calculating the threshold for residues to be considered stabilising
mean_stability <- mean(PDB_mean_per_res$sliding_window, na.rm = TRUE)
SD_stability <- sd(PDB_mean_per_res$sliding_window, na.rm = TRUE)
stabilising_threshold <- mean_stability - SD_stability

# Testing setting 0 as stabilising threshold as -ve dG are stabilising
#stabilising_threshold <- 0

### Visualization ###

# Density Plot
ggplot(data = PDB_mean_per_res, aes(x = sliding_window)) +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = mean_stability, colour = "black", linewidth = 1, alpha = 0.9) +
  geom_vline(xintercept = stabilising_threshold, colour = "red", linewidth = 1, linetype = "dashed", alpha = 0.9) +
  ylab("Density") +
  xlab("Sliding Window Averaged\nMean \u0394G per Residue for Each PDB (kcal.mol-1)") +
  scale_x_continuous(breaks = seq(-10,10, by = 1)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none")

# Saving density plot 
#ggsave(plot = last_plot(), file = "Output/stable_regions/density_plot_with_threshold.png", height = 4, width = 6)

# Line Plot
ggplot(data = PDB_mean_per_res, aes(x=Pos)) +
  geom_line(aes(y = sliding_window, group = interaction(pdb_id, cumsum(c(0, diff(Pos) != 1)))), 
            linewidth = 0.75, alpha = 0.25) +
  geom_hline(yintercept = mean_stability, colour = "black", linewidth = 0.75, alpha = 0.9) +
  geom_hline(yintercept = stabilising_threshold, colour = "red", linewidth = 0.75, linetype = "dashed", alpha = 0.9) +
  ggtitle(paste0("Sliding Window of ", window_size)) +
  ylab("Mean \u0394G per residue (kcal.mol-1)") +
  xlab("Residue Number") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.25, linetype = "dashed"),
        plot.title = element_text(size = 14, colour = "black", face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none")

# Saving density plot
#ggsave(plot = last_plot(), file = "Output/stable_regions/line_plot_with_threshold.png", height = 4, width = 8)

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


# Importing the mean delta G per residue for each PDB 
PDB_mean_per_res <- read.csv("Output/thermodynamics/foldx_stability_filtered_normalised.csv")

# Playing with different sliding windows

window_size <- 3
#window_size <- 5

PDB_mean_per_res <- PDB_mean_per_res %>%
  group_by(pdb_id) %>%
  mutate(sliding_window = rollapply(mean_energy, width = window_size, FUN = mean, align = "center", fill = NA)) %>% 
  ungroup()

# Finding local minima (troughs)
tmp2 <- PDB_mean_per_res %>% group_by(Pos) %>% summarise(mean_Pos_sliding_window = mean(sliding_window, na.rm = TRUE))
tmp3 <- na.omit(tmp2)


mean_stability <- mean(PDB_mean_per_res$sliding_window, na.rm = TRUE)

# Get the turning points object
tp <- pastecs::turnpoints(na.omit(tmp2$mean_Pos_sliding_window), calc.proba = TRUE)

# Creating a turnpoint df from tp
info <- extract(tp, peak = 1, pit = -1)
type <- info[info != 0]

tp_df <- data.frame(Pos = tp$tppos + (min(tmp3$Pos, na.rm = TRUE) - 1), # tppos counts each data point, since my data doesnt start at 1 I need to account for it
                    info = as.character(type),
                    probability = tp$proba)

# Define a probability threshold
threshold <- 0.15  # for instance, if you're interested in low probability values

# Filter the turning points object based on probability
tp_df <- tp_df %>% filter(probability < threshold)

# Create a data frame for plotting
tmp3 <- merge(tmp3, tp_df, by = "Pos", all = TRUE)
tmp3$info[is.na(tmp3$info)] <- "0"

################################################################################


# Visualising peaks
ggplot(tmp3, aes(x = Pos, y = mean_Pos_sliding_window)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_point(aes(color = info), size = 2) +
  scale_color_manual(values = c("1" = "blue", "0" = "black", "-1" = "red"),
                     label= c("Trough", "NA", "Peak")) +
  scale_x_continuous(breaks = seq(0,1000,10)) +
  geom_vline(xintercept = tmp3$Pos[tmp3$info == "1"], color = "blue", linetype = "dashed", size = 0.25) +  # Peaks
  geom_vline(xintercept = tmp3$Pos[tmp3$info == "-1"], color = "red", linetype = "dashed", size = 0.25) +  # Pits
  geom_hline(yintercept = mean_stability, linetype = "solid", size = 1) +  # Mean line
  labs(title = "Turning Points", x = "Residue Number", y = "Mean Pos Sliding Window", colour = "Type") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_line(colour = "grey", linewidth = 0, linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0, linetype = "dashed"),
        plot.title = element_text(size = 14, colour = "black", face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"))

# Saving Turning Points Lineplot
ggsave(plot = last_plot(), file = "Output/stable_regions/turning_points.png", height = 4, width = 8, bg = "white")

#####################################################################################################

### Selecting Stable Regions ###

# Selecting peaks above the mean line
#significant_peaks <- tmp3$Pos[tmp3$info == 1 & tmp3$mean_Pos_sliding_window > mean_stability]
significant_peaks <- tmp3$Pos[tmp3$info == 1]

# Initialize variables
stabilisers <- c()
region_list <- c()
region <- 1
stabiliser_count <- 0

# Loop through the data to find stabilizers
for (i in seq_along(tmp3$mean_Pos_sliding_window)) {
  
  # Check if the current position is stabilizing
  if (tmp3$mean_Pos_sliding_window[i] < mean_stability & tmp3$info[i] != 1) {
    stabilisers <- c(stabilisers, tmp3$Pos[i])
    region_list <- c(region_list, region)
    stabiliser_count <- stabiliser_count + 1
  }
  
  # If a significant peak is hit, check if the current region has more than 1 stabilizer
  if (tmp3$Pos[i] %in% significant_peaks) {
    
    # Only proceed if there are more than 1 stabilizers in the region
    if (stabiliser_count > 1) {
      region <- region + 1
    } else {
      # Remove stabilisers with less than 2 occurrences
      stabilisers <- stabilisers[region_list != region]
      region_list <- region_list[region_list != region]
    }
    
    # Reset stabiliser count for the next region
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
ggplot(tmp3, aes(x = Pos, y = mean_Pos_sliding_window)) +
  geom_rect(data = stable_regions, aes(xmin = min - 0.5, xmax = max + 0.5, ymin = -Inf, ymax = Inf, fill = region),
            alpha = 0.75, inherit.aes = FALSE)+#, fill = "grey50") +
  geom_line(color = "black", linewidth = 1) +
  geom_point(aes(color = info), size = 2) +
  scale_fill_manual(values = region_colours) +
  #scale_fill_brewer(palette = "Set1") +
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

# Saving Turning Points Lineplot
ggsave(plot = last_plot(), file = "Output/stable_regions/turning_points_and_regions.png", height = 5, width = 8, bg = "white")


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

# Saving Turning Points Lineplot
ggsave(plot = last_plot(), file = "Output/stable_regions/sliding_window_and_regions.png", height = 4, width = 8, bg = "white")


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

# Saving Turning Points Lineplot
ggsave(plot = last_plot(), file = "Output/stable_regions/raw_delta_G_and_regions.png", height = 4, width = 8, bg = "white")


# Saving stable regions
stable_regions_minimal <- stable_regions %>% select("residues", "region")

write.csv(stable_regions_minimal, file = "Output/stable_regions/stable_regions.csv", row.names = FALSE)











