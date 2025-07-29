cat("Running: stable_region_distances.R\n")

# Installing and Loading Packages
required_packages <- c("dplyr", "stringr", "readr", "ggplot2", "ggpubr", "bio3d", "tidyr", "ggtext", "igraph", "dendextend")

not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]   

# Install missing packages
if(length(not_installed)) install.packages(not_installed)    

# Load all required packages
invisible(lapply(required_packages, require, character.only = TRUE))

# Cleaning Environment
rm(list = ls())

################################################################################

# Creating Directories
save_dir <- file.path("Output", "stable_regions")
region_dir <- file.path(save_dir, "region_contacts")
residue_dir <- file.path(save_dir, "residue_contacts")

if(!dir.exists(save_dir)) dir.create(save_dir)
if(!dir.exists(region_dir)) dir.create(region_dir)
if(!dir.exists(residue_dir)) dir.create(residue_dir)

################################################################################

# Locate all pdbs in folder that need to be looped through
folder_path <- file.path("Output", "PDBs", "asymetric_unit")

# Get a list of all .pdb files in the folder
pdb_files <- list.files(path = folder_path, pattern = "\\.pdb$", full.names = TRUE)

################################################################################

# Removing PDBs with a poor mean resolution - identified in the validation scripts

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

################################################################################

# Loading Fibril info
fibril_tmp <- read.csv(file.path("Output", "PDBs", "COM_and_fibril.csv"))
fibril_info <- fibril_tmp %>% select(-c("chain", "x_com", "y_com", "z_com", "Rg")) %>% distinct()

# Loading RMSD Cluster Data
cluster_groups <- read.csv(file.path("Output", "RMSD", "data", "RMSD_cluster_groups.csv"))

fibril_df <- merge(fibril_info, cluster_groups, by = "pdb_id")

# Importing stable region data
stable_regions <- read.csv(file.path("Output", "stable_regions", "stable_regions.csv"))

# Separating residues into start and end
stable_regions <- stable_regions %>%
  rowwise() %>%
  mutate(start = as.numeric(str_split_1(residues, "-")[1]),
         end = as.numeric(str_split_1(residues, "-")[2])) %>%
  ungroup()


# Initialising dataframe
combined_residue_distance <- data.frame()
combined_region_distance <- data.frame()

# Defining functions
stable_region_distance <- function(CA_filtered) {
  
  #######################################################
  ##### Calculating Distance Between Stable Regions #####
  #######################################################
  
  CA_filtered$region1 <- NA
  CA_filtered$region2 <- NA
  
  # Finding CA-CA distance comparisons that occur within stable regions
  for (j in  1:nrow(stable_regions)) {
    
    for (k in 1:nrow(CA_filtered)) {
      
      if (CA_filtered$resno1[k] >= stable_regions$start[j] & 
          CA_filtered$resno1[k] <= stable_regions$end[j]) {
        
        CA_filtered$region1[k] <- stable_regions$region[j] 
        
      }
      
      if (CA_filtered$resno2[k] >= stable_regions$start[j] & 
          CA_filtered$resno2[k] <= stable_regions$end[j]) {
        
        CA_filtered$region2[k] <- stable_regions$region[j] 
        
      }
    }
  }
  
  # Removing Non-stable region comparisons
  stable_region_CA_dist <- na.omit(CA_filtered)
  
  # Removing comparisons between same region on same protofibril
  stable_region_CA_dist <- stable_region_CA_dist %>%
    filter(!(region1 == region2 & which_fibril == "same"))
  
  
  # Find the minimum distance for each pair of regions for each fibril
  min_distances <- stable_region_CA_dist %>%
    #group_by(region1, region2) %>%
    group_by(region1, region2, which_fibril, pdb_id) %>%
    summarize(min_distance = min(distance), .groups = 'drop')
  
  
  # Making group columns as characters
  min_distances$region1 <- as.character(min_distances$region1)
  min_distances$region2 <- as.character(min_distances$region2)
  
  # Adding PDB name
  min_distances$pdb <- selected_pdb
  
  return(min_distances)
}


for (file in pdb_files_filtered) {
  
  # Extracting PDB code from the file path
  tmp = tail(str_split_1(file, "/"), 1)
  tmp = str_split_1(tmp, "\\.")[1]
  selected_pdb = str_split_1(tmp, "_")[1]
  
  # Printing PDB Name
  cat("\nAnalysing: ", selected_pdb, "\n")
  
  ##- Read PDB file
  pdb <- read.pdb(file)
  
  # Extracting atom info into a dataframe
  df <- as.data.frame(pdb$atom)
  
  # Since b-factor has been edited to contain the fibril number - changing the name accordingly
  names(df)[names(df) == "b"] <- "fibril"
  
  # adding in the current pdb
  df$PDB <- selected_pdb

  #################################################
  ##### Calculating Distance Between C-alphas #####
  #################################################
  
  # selecting only CA
  df <- df %>% filter(df$elety == "CA")
  
  # Selecting only the relevant columns 
  df <- df %>% select(c("PDB", "fibril", "resno", "x", "y", "z"))
  
  # Adding in pdb_id to handle fibrils with multiple polymorphs
  df <- merge(df, fibril_df, by = c("PDB", "fibril"))
  
  # Counting the number of polymorphs
  polymorph_count <- as.numeric(max(df$polymorph))
  
  
  # Calculate the pairwise distances between CA's using the dist function
  CA_distances <- as.matrix(dist(df[, c("x", "y", "z")]))
  
  # Convert the distance matrix to a long format dataframe
  CA_distance_df <- data.frame(
    resno1 = rep(df$resno, each = nrow(df)),
    fibril1 = rep(df$fibril, each = nrow(df)),
    polymorph1 = rep(df$polymorph, each = nrow(df)),
    resno2 = rep(df$resno, times = nrow(df)),
    fibril2 = rep(df$fibril, times = nrow(df)),
    polymorph2 = rep(df$polymorph, times = nrow(df)),
    distance = as.vector(CA_distances))
  
  # Removing comparisons between same region on same protofibril
  CA_distance_df <- CA_distance_df %>%
    filter(!(resno1 == resno2 & fibril1 == fibril2))
  
  # Handling neighboring residues on the same monomer (removing comparisons with up to 3 neighbouring residues)
  CA_distance_df <- CA_distance_df %>%
    filter(!(abs(resno1 - resno2) <= 3 & fibril1 == fibril2))
  
  # Assigning if the distance calculation is inter- or intra-fibril
  CA_distance_df <- CA_distance_df %>%
    mutate(which_fibril = ifelse(fibril1 == fibril2, "same", "different"))
  
  
  ### If only one unique value in pdb_id (i.e. only one polymorph) ###
  if (polymorph_count == 1) {
    
    # Find the minimum distance for each pair of residues for each fibril (to remove duplicate intra-residue comparisons)
    min_distances <- CA_distance_df %>%
      group_by(resno1, resno2) %>%
      summarize(min_distance = min(distance), .groups = 'drop')
    
    # Join the summarized result back to the original dataframe to retain all columns
    CA_filtered <- CA_distance_df %>%
      inner_join(min_distances, by = c("resno1", "resno2")) %>%
      filter(distance == min_distance) %>%
      select(-c("min_distance", "fibril1", "fibril2"))
    
    if(nrow(CA_filtered > 0)) {
      # Adding pdb_id
      current_pdb_info <- fibril_df %>% filter(PDB == selected_pdb)
      CA_filtered$pdb_id <- current_pdb_info$pdb_id[1]
      
    }
    
    # Merging with combined dataframe to store all pdb distances
    combined_residue_distance <- rbind(combined_residue_distance, CA_filtered)
    
    # Calculating stable region distances
    min_region_distances <- stable_region_distance(CA_filtered)
    
    # Merging current min_unique_data with combined_region_distance
    combined_region_distance <- rbind(combined_region_distance, min_region_distances)
    
  
  } else {
    
    # If there is more than one pdb_id (i.e. more than one polymorph)
    
    # Finding the length of remaining unique polymorphs as some may be removed by Q-score filtering
    remaining_polymorphs <- unique(CA_distance_df$polymorph1)
    
    #for (i in 1:polymorph_count) {
    for (i in remaining_polymorphs) {
        
      # Filter rows where fibril1 or fibril2 contains the polymorph
      CA_distance_polymorph <- CA_distance_df %>%
                                filter(polymorph1 == i)
      
      # Find the minimum distance for each pair of residues for each fibril (to remove duplicate intra-residue comparisons)
      min_distances <- CA_distance_polymorph %>%
        group_by(resno1, resno2) %>%
        summarize(min_distance = min(distance), .groups = 'drop')
      
      # Join the summarized result back to the original dataframe to retain all columns
      CA_filtered <- CA_distance_polymorph %>%
        inner_join(min_distances, by = c("resno1", "resno2")) %>%
        filter(distance == min_distance) %>%
        select(-c("min_distance", "fibril1", "fibril2"))
      
        if(nrow(CA_filtered > 0)) {
          # Getting the polymorph information for the current fibril
          current_pdb_info <- fibril_df %>% filter(PDB == selected_pdb)
          
          # Handling asym unit containing multiple copies of the same polymorph
          current_pdb_info <- current_pdb_info %>%
                                distinct(pdb_id, polymorph, .keep_all = TRUE)
          
          # Dropping fibril column as it may be incorrect for asym_units with multiple copies of the same polymorph
          current_pdb_info <- current_pdb_info %>% select(-c("fibril"))
          
          # Adding pdb_id
          CA_filtered$pdb_id <- current_pdb_info$pdb_id[current_pdb_info$polymorph == i]
          
          }
      
        # Merging with combined dataframe to store all pdb distances
        combined_residue_distance <- rbind(combined_residue_distance, CA_filtered)
        
        # Calculating stable region distances
        min_region_distances <- stable_region_distance(CA_filtered)
        
        # Merging current min_unique_data with combined_region_distance
        combined_region_distance <- rbind(combined_region_distance, min_region_distances)
      
      }
  }
  
}

# Saving combined_residue_distances
write_csv(combined_residue_distance, file.path(residue_dir, "residue_distances.csv"))


# Saving combined_region_distances
write_csv(combined_region_distance, file.path(region_dir, "region_distances.csv"))


#################################################################################################
#################################################################################################
#################################################################################################

##### Creating Network Plots for Stable Region Distances for Each RMSD Cluster Group #####

#################################################################################################
#################################################################################################
#################################################################################################

# Creating a single comparison column
combined_region_distance <- combined_region_distance %>%
  mutate(comparison = paste(region1, "-", region2))

# Merging stable region distances and cluster group dataframes
stable_dist_by_cluster <- merge(combined_region_distance, cluster_groups, by = "pdb_id")

# Saving data
write.csv(stable_dist_by_cluster, file = file.path(region_dir, "region_distances_by_cluster_group.csv"), 
          row.names = FALSE, quote = TRUE)


# Calculating the percentage of structures per group with each contact

# Keeping contacts < 10.8A
contact_df <- stable_dist_by_cluster %>% filter(min_distance < 10.8)

# Counting the number of each contact per group
contact_counts <- contact_df %>%
  count(group, comparison) %>%
  arrange(group)

# Counting the number of unique PDBs in each group
pdbs_per_group <- stable_dist_by_cluster %>%
  group_by(group) %>%
  summarise(unique_pdb_count = n_distinct(pdb_id))

percent_df <- merge(contact_counts, pdbs_per_group, by = "group")

percent_df$freq <- (percent_df$n / percent_df$unique_pdb_count) * 100

percent_df$comparison <- gsub("-", "--", percent_df$comparison)

write.csv(percent_df, file = file.path(region_dir, "region_contacts_percentage_per_group.csv"), 
          row.names = FALSE)


#####################
### Network Plots ###
#####################

# Function to calculate node angles with one node fixed at 3 o'clock (igraph always places a node a 3 o'clock)
calculate_node_angles_fixed_at_3oclock <- function(num_nodes) {
  if (num_nodes <= 0) {
    stop("Number of nodes must be greater than 0")
  }
  
  angle_between_nodes <- 360 / num_nodes
  # Adjust starting angle to fix one node at 3 o'clock (0 degrees)
  starting_angle <- 0
  node_angles <- starting_angle + seq(0, by = angle_between_nodes, length.out = num_nodes)
  node_angles <- node_angles %% 360
  
  return(node_angles)
}

# Function to calculate perpendicular angles to the tangent for each node
calculate_perpendicular_angles <- function(num_nodes) {
  angles <- calculate_node_angles_fixed_at_3oclock(num_nodes)
  
  # Calculate perpendicular angles (90 degrees offset from node angle)
  perpendicular_angles <- (angles + 90) %% 360
  
  # Create a data frame for the results
  result_df <- data.frame(
    Node_Angle = angles,
    Perpendicular_Angle = perpendicular_angles
  )
  
  return(result_df$Perpendicular_Angle)
}

# Calculate the angle to place the self loops so the face outside the circle
num_nodes <- as.numeric(nrow(stable_regions)) 
angles <- calculate_perpendicular_angles(num_nodes)


# Setting distance threshold
distance_threshold = 10.8

network_data <- stable_dist_by_cluster %>% filter(min_distance <= distance_threshold)

# finding all the unique cluster groups
groups <- unique(cluster_groups$group)

for (current_group in groups) {

  # Filtering stable_dist_by_cluster to select only one cluster group at a time
  network_data <- stable_dist_by_cluster %>% 
                    filter(min_distance <= distance_threshold & group == current_group)
  
  if (nrow(network_data) > 0) {
    
    # Getting vertex names and sorting in descending order
    vertex_names <- sort(unique(c(stable_dist_by_cluster$region1, stable_dist_by_cluster$region2)), decreasing = TRUE)
    
    # Create a complete graph with all nodes
    all_nodes <- vertex_names
    all_nodes <- as.character(all_nodes)
    complete_graph <- graph.full(length(all_nodes))
    V(complete_graph)$name <- all_nodes
    
    # Remove all edges from the complete graph
    complete_graph <- delete_edges(complete_graph, E(complete_graph))
    
    
    # Summarize edge weights by counting occurrences
    network_data <- network_data %>%
      group_by(region1, region2, which_fibril) %>%
      summarise(weight = n(), .groups = 'drop')  # Count number of occurrences
    
    # Create the graph using weighted edges
    filtered_graph <- graph_from_data_frame(network_data, directed = FALSE)
    
    # Assign weights to the edges
    E(filtered_graph)$weight <- network_data$weight
    
    # Use weight to determine edge width (scale as needed)
    E(filtered_graph)$width <- 1 + E(filtered_graph)$weight * 0.3  # Adjust multiplier as needed
    
    
    # Merge the complete graph with the filtered graph
    graph <- union(complete_graph, filtered_graph)
    
    # Identify self-loops and color them red
    self_loops <- which_loop(graph)
    
    # Match the sorted edges with network_data
    network_data <- network_data[order(network_data$region2, network_data$region1), ]
    
    # Map which_fibril to colours
    edge_colours <- ifelse(network_data$which_fibril == "same", "grey30", "red")
    
    # Assign these colors to the edges
    E(graph)$color <- edge_colours

    # Calculate the degree of each vertex
    vertex_degrees <- degree(graph)
    
    # Set vertex size based on degree (scaling it for better visualization)
    vertex_sizes <- vertex_degrees * 3
    
    # Extending colours as necessary
    # Define your region size
    region <- nrow(stable_regions)
    
    # Define the original named color vector
    colors <- c("1" = "red", "2" = "blue", "3" = "gold2", "4" = "purple", 
                "5" = "darkorange2", "6" = "cyan", "7" = "green", "8" = "pink",
                "9" = "chocolate4", "10" = "darkgreen", "11" = "magenta", "12" = "deepskyblue",
                "13" = "darkviolet", "14" = "coral", "15" = "darkslategray")
    
    # Generate extended keys for the new vector
    new_keys <- as.character(seq_len(region))
    
    # Extend the colors to match the region size
    node_colors <- rep(colors, length.out = region)
    names(node_colors) <- new_keys
    
    # Assign colors to the vertices
    V(graph)$color <- node_colors[V(graph)$name]
    
    #########################################################
    #########################################################
    #########################################################
    
    # Adding transparency to the nodes
    # Function to convert named colors to RGBA
    add_alpha <- function(col, alpha = 0.5) {
      # Convert color names to RGB and then apply alpha
      rgb_vals <- col2rgb(col) / 255
      apply(rgb_vals, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
    }
    
    # Apply transparency to your node colors
    V(graph)$color <- add_alpha(node_colors[V(graph)$name], alpha = 0.25)
    
    # Changing the border colour
    V(graph)$frame.color <- node_colors[V(graph)$name]
    
    #########################################################
    #########################################################
    #########################################################
    
    # Set border width for each node
    vertex_border_width <- 10  # Adjust this value to control the border width
    
    
    # Define a function to draw a circle
    draw_circle <- function(center, radius, npoints = 100) {
      theta <- seq(0, 2 * pi, length.out = npoints)
      x <- radius * cos(theta) + center[1]
      y <- radius * sin(theta) + center[2]
      return(data.frame(x, y))
    }
    
    # Calculate the center and radius for the circle
    layout <- layout_in_circle(graph)
    center <- colMeans(layout)
    radius <- max(sqrt(rowSums((layout - center)^2)))
    
    # Generate points for the circle
    circle_points <- draw_circle(center, radius)
    
    # Calculate the plot limits to ensure nodes are not clipped
    node_radius <- 30  # Same as vertex.size
    
    border_margin <- vertex_border_width
    total_margin <- node_radius + border_margin
    
    # Calculate plot limits based on layout and margins
    xlim <- range(layout[,1]) + c(-total_margin, total_margin) * 0.015  
    ylim <- range(layout[,2]) + c(-total_margin, total_margin) * 0.015  
    
    # Adjusting the angle of self-loops
    edgeloopAngles <- numeric(ecount(graph))
    b <- 1
    m <- 0
    
    # Create a count vector with length equal to the number of vertices
    node_counts <- unique(stable_regions$region)
    
    angles <- 2 * pi * (seq(num_nodes - 1, 0, by = -1) / num_nodes)
    
    # Loop through edges
    for (edge_id in 1:ecount(graph)) {
      ends <- ends(graph, edge_id)
      v <- which(V(graph)$name == ends[1])
      u <- which(V(graph)$name == ends[2])
      
      if (v == u) {
        m <- m + 1
        
        if (v <= length(angles)) {
          angle <- angles[v]
          edgeloopAngles[b] <- angle + (0.05 * node_counts[v])
          node_counts[v] <- node_counts[v] + 1
        } else {
          stop("Vertex number exceeds the predefined angles list")
        }
        
      } else {
        edgeloopAngles[b] <- 0  # Placeholder for non-self-loops
      }
      
      b <- b + 1
    }
    
    # Open a PNG graphics device
    png(filename = file.path(region_dir, paste0("group_", current_group, "_network.png")),
        width = 10, height = 8, units = "in", res = 300, bg = "transparent")  

    # Plot network graph
    plot(graph, 
         layout = layout,
         edge.width = E(graph)$width,
         #edge.curved = edge_curvatures,
         edge.loop.angle = edgeloopAngles,
         vertex.size = node_radius,
         vertex.label.cex = 2,
         vertex.label.font = 2,
         vertex.label.color = "black",
         vertex.frame.width = vertex_border_width, 
         rescale = FALSE,
         xlim = xlim, ylim = ylim)
    
    mtext(paste0("Group ", current_group), line = -2, cex = 3, font = 2)
    
    # Close the PNG graphics device
    dev.off()
    
  }
}

####################################################################################
####################################################################################
####################################################################################

##### Creating Network Plots for Residue Distances for Each RMSD Cluster Group #####

####################################################################################
####################################################################################
####################################################################################

### Adding RMSD cluster groups ###

# Removing _chainID to just use PDB codes 
cluster_groups <- cluster_groups %>%
  mutate(pdb = sapply(str_split(pdb_id, "_"), function(x) x[1]))

# Merging stable region distances and cluster group dataframes
residue_distance_by_cluster <- merge(combined_residue_distance, cluster_groups, by = "pdb_id")

# Get unique nodes in the entire range (even those without edges)
all_nodes <- seq(min(residue_distance_by_cluster$resno1), max(residue_distance_by_cluster$resno2), 1)

# Getting a list of the groups
groups <- unique(residue_distance_by_cluster$group)

for (current_group in groups) {
  
  tmp <- residue_distance_by_cluster %>% filter(distance <= distance_threshold, group == current_group)
  
  if (nrow(tmp) > 0) {
    
    # Count occurrences of each edge (weight calculation)
    tmp <- tmp %>%
      group_by(resno1, resno2, which_fibril) %>%
      summarise(weight = n(), .groups = 'drop')

    # Create an empty graph with all nodes
    graph <- make_empty_graph(n = length(all_nodes), directed = FALSE)
    V(graph)$name <- all_nodes
    
    # Add edges from the filtered dataframe
    for(i in 1:nrow(tmp)) {
      graph <- add_edges(graph, c(as.character(tmp$resno1[i]), as.character(tmp$resno2[i])))
    }
    
    # Add edges with weights
    E(graph)$weight <- tmp$weight  # Assign weights
    
    # Use weight to determine edge width (scale as needed)
    E(graph)$width <- 0.75 + E(graph)$weight * 0.05  # Adjust multiplier as needed
    
    
    # Map which_fibril to colors
    edge_colours <- ifelse(tmp$which_fibril == "same", "grey30", "red")
    
    # Calculate curvature based on the numeric difference between nodes
    edge_curvature <- ifelse(abs(tmp$resno1 - tmp$resno2) < 20, -0.5, 0)
    edge_curvature[abs(tmp$resno1 - tmp$resno2) < 10] <- -1
    
    
    # Create a circular layout and order it
    layout_sorted <- layout_in_circle(graph)
    
    # Order nodes numerically in an anticlockwise fashion - to match region plot
    node_order <- order(-as.numeric(V(graph)$name))  # Negative ensures decreasing order
    layout_sorted <- layout_sorted[node_order, ]
    
    
    ### Colour nodes to match with stable regions ###
    
    # Separating residues into start and end
    stable_regions <- stable_regions %>%
                        rowwise() %>%
                          mutate(start = as.numeric(str_split_1(residues, "-")[1]),
                                 end = as.numeric(str_split_1(residues, "-")[2])) %>%
                            ungroup()
    
    # Initialize node colors with a default color
    node_colours <- rep("white", length(all_nodes))
    
    region_colours <- c("red", "blue", "gold2", "purple", "darkorange2", "cyan", "green", "pink", "chocolate4", 
                        "darkgreen", "magenta", "deepskyblue", "darkviolet", "coral", "darkslategray")
    
    j = 0
    
    # Assign colors to nodes based on stable regions
    for (i in 1:nrow(stable_regions)) {
      region_start <- stable_regions$start[i]
      region_end <- stable_regions$end[i]
      # Cycle through colours
      region_colour <- region_colours[j %% length(region_colours) + 1]
      
      j = j + 1
      
      # Set the color for nodes within the stable region range
      node_colours[all_nodes >= region_start & all_nodes <= region_end] <- region_colour
    }
    
    
    # Plot network graph
    plot(graph,
         layout = layout_sorted,
         edge.color = edge_colours,
         edge.width = E(graph)$width,
         edge.curved = edge_curvature,
         vertex.size = 5,
         vertex.label.cex = 1,
         vertex.label.font = 2,
         vertex.label.color = "black",
         vertex.frame.width = 2,
         vertex.color = node_colours)
    
    
    # Open a PNG graphics device
    png(filename = file.path(residue_dir, paste0("group_", current_group, "_residue_distances.png")),
        width = 15, height = 15, units = "in", res = 300, bg = "transparent")  # Adjust width, height, and resolution as needed
    
      plot(graph,
           layout = layout_sorted,
           edge.color = edge_colours,
           edge.width = E(graph)$width,
           edge.curved = edge_curvature,
           vertex.size = 5,
           vertex.label.cex = 1,
           vertex.label.font = 2,
           vertex.label.color = "black",
           vertex.frame.width = 2,
           vertex.color = node_colours)
      
      mtext(paste0("RMSD Cluster Group ", current_group), line = 0, cex = 3, font = 2)
    
    # Close the PNG graphics device
    dev.off()
    
  }
}
    