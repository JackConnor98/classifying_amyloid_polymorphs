cat("Running: plot_ordered_residues.R\n")

# Installing and Loading Packages
required_packages <- c("dplyr", "stringr", "tidyr", "readr", "ggplot2", "purrr", "bio3d")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]   

# Install missing packages
if(length(not_installed)) install.packages(not_installed)    

# Load all required packages
invisible(lapply(required_packages, require, character.only = TRUE))

# Cleaning Environment
rm(list = ls())

############################################################################################

# Read in PDB metadata
data <- read_delim("Output/selected_pdbs_metadata.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

# Removing NMR structures
data <- data %>% filter(Method != "ssNMR")

# Selecting only PDB and Residues columns
residues <- subset(data, select = c("PDB ID", "Residues Ordered"))

# Renaming variables for ease of use
colnames(residues) <- c("pdb", "ordered_residues")

# Splitting Residues Ordered by "," as these indicate gaps within the chain
residues$ordered_residues <- str_split(residues$ordered_residues, ",")

# Assigning a unique group ID to each stretch of continuous residues in the pdb
# If for one PDB residues 1-10 and 20-60 are solved they will be given different groups
# This is so that the plot can draw lines to continuous regions and leave spaces where residues aren't ordered
residues <- residues %>%
              mutate(dash_count = str_count(ordered_residues, "-"),
                     group = map2(dash_count, 1:dash_count, seq), 
                     group = str_split(group, ":")) %>%
              select(pdb, ordered_residues, group)

# Unnesting the split character vectors into seperate rows and removing c("") symbols left over
residues <- residues %>% 
              unnest(cols = c(ordered_residues, group)) %>% 
                mutate(ordered_residues = str_replace_all(ordered_residues, 
                                                             c("c" = "", "\"" = "", "\\(" = "", "\\)" = "")))

# Splitting each continuous residue group to give a min and max range
# e.g. 20-96 gets split into two rows; 20 and 96 
# On the plot a line will be draw connecting these two dots
residues$ordered_residues <- str_split(residues$ordered_residues, "-")

# Unnesting into seperate rows
df_unnested <- residues %>% unnest_longer(ordered_residues)

# Setting ordered_residues to be numeric
df_unnested$ordered_residues <- as.numeric(df_unnested$ordered_residues)


# Adding in any local PDB files

pdb_files <- list.files(pattern = "\\.pdb$")

if (length(pdb_files) > 0) {

  for (i in 1:length(pdb_files)) {
    
    current_pdb <- read.pdb(paste(pdb_files[i]))
    
    residue_numbers <- current_pdb$atom$resno
    
    min_resno <- min(residue_numbers)
    
    max_resno <- max(residue_numbers)
    
    pdb_name <- gsub(".pdb", "", pdb_files[i])
    
    tmp <- data.frame(pdb = c(pdb_name, pdb_name),
                      ordered_residues = c(min_resno, max_resno), 
                       group = c(1, 1)  
                      )
    
    df_unnested <- rbind(df_unnested, tmp)
    
  }
}

################################################################################

# Getting number of unique PDBs
num_pdbs <- as.numeric(nrow(data))

# Plotting
p <- ggplot(df_unnested, aes(x = pdb, y = ordered_residues)) 
  
  if (any(grepl("synuclein", data$Protein))) {
    
    q <- p + 
      
      # N-term
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 60, 
               alpha = 0.5, fill =  "darkblue", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 7, ymin = 13, ymax = 23, 
               alpha = 0.5, fill =  "darkblue", colour = "black") +
      
      annotate(geom = "text", label = "N-Term", size = 5, 
               x = num_pdbs + 4, y = 18) +
      
      # NAC
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 61, ymax = 95, 
               alpha = 0.5, fill =  "red", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 7, ymin = 73, ymax = 83, 
               alpha = 0.5, fill =  "red", colour = "black") +
      
      annotate(geom = "text", label = "NAC", size = 5, 
               x = num_pdbs + 4, y = 78) +
      
      
      # C-term
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 96, ymax = 140, 
               alpha = 0.5, fill =  "green", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 7, ymin = 113, ymax = 123, 
               alpha = 0.5, fill =  "green", colour = "black") +
      
      annotate(geom = "text", label = "C-Term", size = 5, 
               x = num_pdbs + 4, y = 118) 
      
      
      # P1
      #annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 36, ymax = 42, 
      #         alpha = 0.5, fill =  "orange", colour = "black") +
      
      #annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 7, ymin = 34, ymax = 44, 
      #         alpha = 0.5, fill =  "orange", colour = "black") +
      
      #annotate(geom = "text", label = "P1", size = 5, 
      #         x = num_pdbs + 4, y = 39) +
      
      
      # P2
      #annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 45, ymax = 57, 
      #         alpha = 0.5, fill =  "green", colour = "black") +
      
      #annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 7, ymin = 46, ymax = 56, 
      #         alpha = 0.5, fill =  "green", colour = "black") +
      
      #annotate(geom = "text", label = "P2", size = 5, 
      #         x = num_pdbs + 4, y = 51) 
    
    
  } else if (any(grepl("Amyloid-", data$Protein))) {
    
    # Reference - https://doi.org/10.1021/jp210019h
    
    q <- p + 
      
      # N-Term
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 16, 
               alpha = 0.5, fill =  "darkblue", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 0.75, xmax = num_pdbs + 2.75, ymin = 6, ymax = 10, 
               alpha = 0.5, fill =  "darkblue", colour = "black") +
      
      annotate(geom = "text", label = "N-Term", size = 4, 
               x = num_pdbs + 1.75, y = 8) +
      
      
      # Central Hydrophobic Core
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 17, ymax = 21, 
               alpha = 0.5, fill =  "black", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 0.75, xmax = num_pdbs + 2.75, ymin = 17, ymax = 21, 
               alpha = 0.5, fill =  "black", colour = "black") +
      
      annotate(geom = "text", label = "Hydrophobic\nCore", size = 4, 
               x = num_pdbs + 1.75, y = 19) +
      
      # Turn Region
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 24, ymax = 27, 
               alpha = 0.5, fill =  "red", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 0.75, xmax = num_pdbs + 2.75, ymin = 23.5, ymax = 27.5, 
               alpha = 0.5, fill =  "red", colour = "black") +
      
      annotate(geom = "text", label = "Turn\nRegion", size = 4, 
               x = num_pdbs + 1.75, y = 25.5) +
      
      
      # Second Hydrophobic Region
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 28, ymax = 35, 
               alpha = 0.5, fill =  "orange", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 0.75, xmax = num_pdbs + 2.75, ymin = 29.5, ymax = 33.5, 
               alpha = 0.5, fill =  "orange", colour = "black") +
      
      annotate(geom = "text", label = "Hydrophobic\nRegion", size = 4, 
               x = num_pdbs + 1.75, y = 31.5) +
      
      
      # C-Term
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 36, ymax = 42, 
               alpha = 0.5, fill =  "green", colour = "black") +
    
    annotate(geom = "rect", xmin = num_pdbs + 0.75, xmax = num_pdbs + 2.75, ymin = 37, ymax = 41, 
             alpha = 0.5, fill =  "green", colour = "black") +
      
      annotate(geom = "text", label = "C-Term", size = 4, 
               x = num_pdbs + 1.75, y = 39) 
    
    
      
  } else if (any(grepl("Tau", data$Protein))) {
    
    # Reference - https://doi.org/10.1021/jp210019h
    
    q <- p + 
      
      # R1
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 243, ymax = 273, 
               alpha = 0.5, fill =  "red", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 5, ymin = 253, ymax = 263, 
               alpha = 0.5, fill =  "red", colour = "black") +
      
      annotate(geom = "text", label = "R1", size = 5, 
               x = num_pdbs + 3, y = 258) +
      
      # R2
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 274, ymax = 304, 
               alpha = 0.5, fill =  "green", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 5, ymin = 284, ymax = 294, 
               alpha = 0.5, fill =  "green", colour = "black") +
      
      annotate(geom = "text", label = "R2", size = 5, 
               x = num_pdbs + 3, y = 289) +
      
      
      # R3
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 305, ymax = 335, 
               alpha = 0.5, fill =  "orange", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 5, ymin = 315, ymax = 325, 
               alpha = 0.5, fill =  "orange", colour = "black") +
      
      annotate(geom = "text", label = "R3", size = 5, 
               x = num_pdbs + 3, y = 320) +
      
      
      # R4
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 336, ymax = 367, 
               alpha = 0.5, fill =  "darkblue", colour = "black") +
      
      annotate(geom = "rect", xmin = num_pdbs + 1, xmax = num_pdbs + 5, ymin = 346.5, ymax = 356.5, 
               alpha = 0.5, fill =  "darkblue", colour = "black") +
      
      annotate(geom = "text", label = "R4", size = 5, 
               x = num_pdbs + 3, y = 351.5) 
      
  } else { 
    q <- p 
  }
  
r <- q + geom_point(size = 3) +
  geom_line(aes(group = pdb), linewidth = 1, linetype = "dotted", colour = "black") +
  geom_line(aes(group = interaction(pdb, group)), linewidth = 2) +
  #ggtitle("Ordered Residues") +
  ylab("Residue Position") +
  xlab("PDB") +
  scale_x_discrete() +
  scale_y_continuous(breaks = seq(0,9999,10)) +
  
  coord_cartesian(xlim = c(1, num_pdbs), clip = "off") +
  
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text.y = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black", face = "bold", angle = 90, hjust=1, vjust = 0.4),
        axis.title=element_text(size = 25, face = "bold"),
        axis.line = element_line(color="black", linewidth = 1.5),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        
        plot.margin = unit(c(1,6,1,1), "lines")) # widens the right margin

plot(r)

# Saving Plot
ggsave(plot = last_plot(), file = "Output/Ordered_Residues.png", height = 8, width = 20)

