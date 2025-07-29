cat("Running: Rg_plotting.R\n")

############################################################################################

# Installing and Loading Packages
required_packages <- c("dplyr", "stringr", "tidyr", "readr", "ggplot2")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]   

# Install missing packages
if(length(not_installed)) install.packages(not_installed)    

# Load all required packages
invisible(lapply(required_packages, require, character.only = TRUE))

# Cleaning Environment
rm(list = ls())

############################################################################################

# Load Data
df <- read_csv(file.path("Output", "PDBs", "COM_and_fibril.csv"))


# Mean Rg for each fibril is calculated
fibril_mean_df <- df %>%
                  group_by(PDB, fibril) %>%
                  summarise(fibril_mean_rg = mean(Rg))


diff_df <- fibril_mean_df %>%
            group_by(PDB) %>%
              mutate(diff_from_fibril1 = (abs(fibril_mean_rg - fibril_mean_rg[fibril == 1]) / fibril_mean_rg) * 100)

# Ordering by difference magnitude
rg_diff_order <- diff_df %>%
  group_by(PDB) %>%
  summarise(max_diff = max(diff_from_fibril1)) %>%
  arrange(max_diff) #arrange(desc(max_diff)) for largest to smallest

# Set the factor levels in that order
diff_df$PDB <- factor(diff_df$PDB, levels = rg_diff_order$PDB)

# Setting fibril number as character for discrete colours
diff_df$fibril <- as.character(diff_df$fibril)

# Plotting
ggplot(data = diff_df, aes(x = PDB, y = diff_from_fibril1, colour = fibril)) +
  geom_point(size = 3, alpha = 0.75,  stroke = 1) +
  geom_hline(yintercept = 5, colour = "black",  linetype = "dashed", linewidth = 1) +
  labs(y = "% difference in the mean Rg\nfor each fibril compared to fibril 1",
       colour = "Fibril Number") +
  theme_bw() +
  theme(axis.title = element_text(size = 16, colour = "black", face = "bold"),
        axis.text = element_text(size = 11, colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.title = element_text(size = 16, colour = "black", face = "bold"),
        legend.text = element_text(size = 14, colour = "black")) 

# Saving plot
ggsave(plot = last_plot(), file = file.path("Output", "PDBs", "rg_difference_plot.png"),
       height = 5, width = 15)

# Saving max Rg difference data
write.csv(rg_diff_order, file = file.path("Output", "PDBs", "maximum_fibril_rg_differences.csv"))
