##################################################################
################## COMPUTE INTERACTIONS ##########################
##################################################################

# it requires a "gap" to be defined, should be required by the function
# Gap is 10 seconds as these are interactions (he maximum gap in tracking before cutting the trajectory or interactions in two different object.)
compute_Interactions <- function(e, start, end, max_time_gap) { # min_cum_duration , link_type, nest_focus, frm_rate
  
  print(paste0("Computing interactions"))
  if (!exists("FRAME_RATE")) {
    stop("FRAME_RATE not defined, can't compute the interactions' duration. Define a frame rate (fps), it will likely be one of these values: 2,4,6,8. ")
  }
  if (!exists("interactions_of_interest")) {
    warning("No fmMatcherInteractionType defined, interactions computed over all of ant shape overlaps combinations. \n If interested in querying interactions only for specific body parts overlaps, list them as follows:\n interactions_of_interest <- c(list(c(\"head\",\"body\")),list(c(\"head\",\"head\")))
")
  }

  # general
  N_DECIMALS <- 3 ## when assigning time in seconds, the number of decimals to preserve when rounding to match between interactions & collisions

  # required packages
  require(FortMyrmidon)

  ### initialise variables
  list_IF_Frames <- list()

  ###############################################################################
  ###### IDENTIFY FRAMES ########################################################
  ###############################################################################

  # Because of various issues raised, including the ones reported here https://github.com/formicidae-tracker/myrmidon/issues/240 ,
  # the analysis is performed not using UNIX_time but using frames. Every new queried file will then show a part including FRAME assignment.
  # Frames are then used for any trajectory cutting later on in the sub-scripts
  IdentifyFrames <- fmQueryIdentifyFrames(e, start = start, end = end, showProgress = FALSE)
  IF_frames <- IdentifyFrames$frames
  # LS: was added because in some cases time was not sorted properly (only a few rows would be off) ####
  IF_frames           <- IF_frames[order(IF_frames$time), ] 
  
  # Assign a frame to each time since start and use it as baseline for all matching and computation
  IF_frames$frame_num <- as.numeric(seq.int(nrow(IF_frames)))

  # assign a time in sec to match annotation time (UNIX hard to use for class mismatch)
  IF_frames$time_sec <- round(as.numeric(IF_frames$time), N_DECIMALS)
  
  
  ### NATH_FLAG
  # deduce FRAME_RATE from IF_frames
  # FRAME_RATE <- rev(c(1:30))[match.closest(x=median(diff(IF_frames$time_sec)),table=rev(1/c(1:30)))]

  # assign frame numbering to annotations
  ### NATH_FLAG: why do you not create a new object that subsets the annotations for the replicate and period you want?
  ### here you have an object with lots of irrelevant rows

  # annotations$frame_start <- match.closest(x = annotations$T_start_sec, table = IF_frames$time_sec, tolerance = 0.05)
  # annotations$frame_stop  <- match.closest(x = annotations$T_stop_sec,  table = IF_frames$time_sec, tolerance = 0.05)
  #
  # Creating a new zeroed-time since the start of the exp by  summing the cumulated differences between each pair of consecutive frames (to account for the time mismatch)
  IF_frames$cum_diff <- c(0, with(IF_frames, diff(time)))
  IF_frames$cum_diff <- cumsum(IF_frames$cum_diff)

  ###############################################################################
  ###### READING TRAJECTORIES ###################################################
  ###############################################################################

  # I believe I don't need this section from BEH_extract_movement_variables_xxxxx.R

  #########################################################################################
  ###### READING AUTOMATIC INTERACTIONS ###################################################
  #########################################################################################

  # Adapted from BEH_Extract_movement_variables_fort081_FUNCTIONS.R

  # Interactions <- fmQueryComputeAntInteractions(e, start, end, maximumGap=gap, singleThreaded=FALSE, reportFullTrajectories = F)

  # Get info on the capules shapes and use the relevant ones in the ant interaction query
  capsules <- e$antShapeTypeNames
  names(capsules) <- as.character(1:length(capsules))

  ALL_CAPS_MATCHERS <- list()
  for (interaction_of_interest in interactions_of_interest) {
    caps_matcher <- c()
    for (caps in interaction_of_interest[c(1:2)]) {
      caps_matcher <- c(caps_matcher, as.numeric(names(capsules)[[which(capsules == caps)]]))
    }
    ALL_CAPS_MATCHERS <- c(ALL_CAPS_MATCHERS, list(caps_matcher))
  }

  # COMPUTE INTERACTIONS, using the same parameters as in Stroeymeyt et al., Science 2018
  # interaction_detection calls merge_interactions.cpp to fix interactions splitting issue
  # interaction_detection automatically performs the analysis on large time chunks (12h)
  # in interaction_detection.R
  Interactions <- interaction_detection(
    e = e,
    start = start,
    end = end,
    max_time_gap = max_time_gap,
    max_distance_moved = (body_lengths$body_length_px / body_lengths$body_length_mm) * MAX_DIST_GAP_MM,
    capsule_matcher = ALL_CAPS_MATCHERS,
    IF_frames = IF_frames
  )

  # ## Add ant_x names and times to the Interactions to convert from FRAME since the start of the experiment, to FRAMES
  # ##creates a ID string for each ant in $interactions
  # Interactions$ant1ID_str            <- paste("ant_",Interactions$ant1,sep="")
  # Interactions$ant2ID_str            <- paste("ant_",Interactions$ant2,sep="")
  # Assign interaction pair
  # Interactions$pair <- paste(Interactions$ant1, Interactions$ant2, sep="_") ## ant 1 is always < ant 2, which makes things easier...

  # ## convert times to different formats
  # Interactions$T_start_UNIX <- as.POSIXct(Interactions$start, format = "%Y-%m-%dT%H:%M:%OSZ", origin = "1970-01-01", tz = "GMT")
  # Interactions$T_stop_UNIX <- as.POSIXct(Interactions$end, format = "%Y-%m-%dT%H:%M:%OSZ", origin = "1970-01-01", tz = "GMT")

  # calc duration (including start frame)
  Interactions$duration <- round((Interactions$end - Interactions$start + 1 / FRAME_RATE), N_DECIMALS)

  #assign time in sec to avoid issues on time management and matching
  Interactions$Starttime <- round(as.numeric(Interactions$start),N_DECIMALS)
  Interactions$Stoptime <- round(as.numeric(Interactions$end),N_DECIMALS)
  
  # round digits
  Interactions$ant1.mean.x          <- round(as.numeric(Interactions$ant1.mean.x),N_DECIMALS)
  Interactions$ant1.mean.y          <- round(as.numeric(Interactions$ant1.mean.y),N_DECIMALS)
  Interactions$ant1.mean.angle      <- round(as.numeric(Interactions$ant1.mean.angle),N_DECIMALS)
  Interactions$ant2.mean.x          <- round(as.numeric(Interactions$ant2.mean.x),N_DECIMALS)
  Interactions$ant2.mean.y          <- round(as.numeric(Interactions$ant2.mean.y),N_DECIMALS)
  Interactions$ant2.mean.angle      <- round(as.numeric(Interactions$ant2.mean.angle),N_DECIMALS)
  
  # check that frame numbering is correct
  # #assign start and end frame number
  # Interactions$frame_start <- match.closest(x = Interactions$T_start_sec, table = IF_frames$time_sec, tolerance = 0.05)
  # Interactions$frame_stop  <- match.closest(x = Interactions$T_stop_sec,  table = IF_frames$time_sec, tolerance = 0.05)
  ### check that frame numbering is correct
  # tail(Interactions$startframe)
  # tail(Interactions$frame_start)

  # Rename variables in a SCIENCE 2018 COMPATIBLE FORMAT
  # EXPECTED FORMAT:: Tag1,Tag2,Startframe,Stopframe,Starttime,Stoptime,Box, Xcoor1,Ycoor1,Angle1,Xcoor2,Ycorr2,Angle2,Direction (types),Detections,
  # DIRECTION → the capsule types matching
  # DETECTION → the rate of frames in which both ants where present
  names(Interactions)[which(names(Interactions) == "ant1")] <- "Tag1" # beware: these are antID but called Tags for compatibility
  names(Interactions)[which(names(Interactions) == "ant2")] <- "Tag2"
  names(Interactions)[which(names(Interactions) == "types")] <- "Direction"
  names(Interactions)[which(names(Interactions) == "detections")] <- "Detections"
  names(Interactions)[which(names(Interactions) == "startframe")] <- "Startframe"
  names(Interactions)[which(names(Interactions) == "endframe")] <- "Stopframe"
  # names(Interactions)[which(names(Interactions) == "T_start_UNIX")] <- "Starttime"
  # names(Interactions)[which(names(Interactions) == "T_stop_UNIX")] <- "Stoptime"
  names(Interactions)[which(names(Interactions) == "ant1.mean.x")] <- "Xcoor1"
  names(Interactions)[which(names(Interactions) == "ant1.mean.y")] <- "Ycoor1"
  names(Interactions)[which(names(Interactions) == "ant1.mean.angle")] <- "Angle1"
  names(Interactions)[which(names(Interactions) == "ant2.mean.x")] <- "Xcoor2" # fixed typo in original naming Xcoor2
  names(Interactions)[which(names(Interactions) == "ant2.mean.y")] <- "Ycoor2"
  names(Interactions)[which(names(Interactions) == "ant2.mean.angle")] <- "Angle2"
  Interactions$Box <- e$spaces[[1]]$name # tracking system
  # remove extras
  Interactions$start          <- NULL
  Interactions$end            <- NULL
  Interactions$space          <- NULL
  Interactions$T_start_UNIX   <- NULL
  Interactions$T_stop_UNIX    <- NULL

  # should I retain the Tag1,Tag2 or is ant1 ant2 enough as an identifier?

  # TRUNCATION
  # truncate interactions over 2 mins
  long_interaction_indices <- which(Interactions$duration > 120)
  Interactions[long_interaction_indices,"duration"] <- 120
  Interactions[long_interaction_indices,"Stoptime"] <-  Interactions[long_interaction_indices,"Starttime"] + 120
  Interactions[long_interaction_indices,"Stopframe"] <-  Interactions[long_interaction_indices,"Startframe"] + round(  ( Interactions[long_interaction_indices,"duration"] - 1/FRAME_RATE )* FRAME_RATE)
  
  
  return(Interactions)
}


#### BEST PRACTICE:
# SEPARATE THE INTERACTION CREATION SCRIPT FROM THE NETWORK SCRIPT AND SAVE ITS OUTPUT.
# IF THE OUTPUT IS ALREADY PRESENT, THEN LOAD THE INTERACTIONS FILE

### THE COMPUTE GRAPH/NETWORK MEASURES SCRIPTS SHOULD MANAGE THE INTERACTION FILE BY
#   DEALING WITH PRE AND POST SEPARATELY (CURRENTLY THEY WORK INSIDE THE 3 H LOOP BUT IT IS PROBABLY
#   RETUNDANT/USELESS BY NOW). reformatting of the output of these functions may be needed.

## The Ant Task assignment  and space use benefit from the time windows, don't mess with them!!!!


##################################################################
################## COMPUTE NETWORK ###############################
##################################################################

compute_G <- function(e, Interactions) {
  print(paste0("Computing networks"))

  # required packages
  require(igraph)

  # Interactions$ant_ID1 <- paste("ant_",Interactions$ant1,sep="") ##creates a ID string for each ant: ant1, ant2,...
  # Interactions$ant_ID2 <- paste("ant_",Interactions$ant2,sep="")
  e.Ants <- e$ants
  Ant_IDs <- NULL
  # Number of ants
  for (ant in 1:length(e.Ants)) {
    Ant_IDs <- c(Ant_IDs, e.Ants[[ant]]$ID)
  }
  Max_antID <- max(Ant_IDs)
  # N_ants <- length(e.Ants)

  # initialise adj-matrix
  adj_matrix <- matrix(0L, nrow = Max_antID, ncol = Max_antID) # np.zeros((N_ants, N_ants))
  rownames(adj_matrix) <- Ant_IDs
  colnames(adj_matrix) <- Ant_IDs
  # Populate network
  for (INT in 1:nrow(Interactions)) {
    ANT1 <- Interactions[INT, "ant1"]
    ANT2 <- Interactions[INT, "ant2"]

    # Populate adjaciency matrix
    # Choose “link as number of interactions”:
    adj_matrix[ANT1, ANT2] <- adj_matrix[ANT1, ANT2] + 1
    # or “link as total duration of interactions”:
    #  adj_matrix[ant_ID1, ant_ID2] = adj_matrix[ant_ID1, ant_ID2] + int$time_stop – int$time_start

    # if link_type == 'length_inter':
    #   # OPT1
    #   # WEIGHTS: cumulative interaction time
    #   adj_mat[i.IDs[0]-1, i.IDs[1]-1] += (TimeToFrame[fm.Time.ToTimestamp(i.End)] - TimeToFrame[fm.Time.ToTimestamp(i.time_start)]) * 1 / frm_rate
    # elif link_type == '#inter':
    #   # OPT2
    #   # WEIGHTS: number of interactions
    #   adj_mat[i.IDs[0]-1, i.IDs[1]-1] += 1
    # else:
    #   raise TypeError('"link_type" not valid')
  }

  adj_matrix <- adj_matrix + t(adj_matrix)

  # interaction filtering (remove weak connections)
  # adj_mat[adj_mat <  min_cum_duration] = 0

  # network build
  G <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  actors <- V(G)
  # # store inverse of weights
  # ADRIANO: look for the function in igraph that works as set_edge_attributes in networkx of python
  # nx.set_edge_attributes(G,
  #                        {(i,j): 1/adj_mat[j,i] if adj_mat[j,i]>0 else 0 for i in range(len(adj_mat)) for j in range(i)},
  #                        'inv_weight')

  #### add a column contaning interaction duration in min
  # Interactions["duration_min"] <- as.numeric(difftime(Interactions$time_stop, Interactions$time_start, units = "mins") + 0.125) ###duration in minutes (one frame = 0.125 second)
  ### add edge weights
  E(G)$weight <- Interactions[, "duration"]
  ### simplify graph (merge all edges involving the same pair of ants into a single one whose weight = sum of these weights)
  G <- simplify(G, remove.multiple = TRUE, remove.loop = TRUE, edge.attr.comb = "sum")
  ################## remove unconnected nodes
  unconnected <- actors[degree(G) == 0]
  G <- G - as.character(unconnected)
  ################## update actor list
  actors <- get.vertex.attribute(G, "name")

  ################# assign vertex attributes
  # subset metadata
  METADATA <- subset(metadata, metadata$REP_treat == REP_TREAT)
  # create matching of vertices with ants (there could be less ants in a 3-hours window than in the full e)

  ############# Assign vertex types: "AntTask_num"
  if ("AntTask" %in% colnames(metadata)) {
    METADATA$AntTask_num <- NA
    METADATA[which(METADATA$AntTask == "nurse"), "AntTask_num"] <- 1
    METADATA[which(METADATA$AntTask == "forager"), "AntTask_num"] <- 2

    # keep only ants corresponding to vertices
    V_METADATA <- METADATA[which(METADATA$antID %in% V(G)$name), ]

    G <- set_vertex_attr(G, name = "AntTask_num", index = V(G), value = V_METADATA$AntTask_num)
  } else {
    warning("No metadata provided, impossible to assign Vertex types to Graph")
  }

  #############  Assign vertex types: "AntTask"
  G <- set_vertex_attr(G, name = "AntTask", index = V(G), value = V_METADATA$AntTask)

  #############  Assign vertex types: "Exposed"
  G <- set_vertex_attr(G, name = "Exposed", index = V(G), value = V_METADATA$Exposed)

  #############  Assign vertex types: "IsQueen"
  G <- set_vertex_attr(G, name = "IsQueen", index = V(G), value = V_METADATA$IsQueen)

  #############  Assign vertex types: "antID"
  G <- set_vertex_attr(G, name = "antID", index = V(G), value = V_METADATA$antID)

  # REMOVE nodes without ANTTASK (missing ants from observation as likely died at early stage of experiment - not appearing in the 48h pre-exposure)
  G <- delete_vertices(G, is.na(V(G)$AntTask_num))

  return(G)
} # compute_G

##################################################################
################## COMPUTE NETWORK PROPERTIES ####################
##################################################################

NetProperties <- function(graph) {
  print(paste0("Computing network properties"))

  # required packages
  require(igraph)
  require(FortMyrmidon)


  summary_collective <- NULL
  summary_individual <- NULL


  ##### COLLECTIVE NETWORK PROPERTIES ######################################
  #### inherited from Stroeymeyt et al. 2018

  ## Assortativity - Task
  # degree of preferential association between workers of the same task group, calculated using Newman’s method
  # if (!is.n(V(graph)$AntTask_num)) {
  task_assortativity <- assortativity_nominal(graph, types = V(graph)$AntTask_num, directed = F) # }else{
  # warning("No  Task Vertex information, impossible to calculate task_assortativity")}

  ## Clustering
  clustering <- mean(transitivity(graph, type = "barrat", weights = E(graph)$weight, isolates = c("NaN")), na.rm = T)
  ## Degree mean and max
  # Degree centrality:   degree of a vertex is its the number of its adjacent edges.
  degrees <- degree(graph, mode = "all")
  degree_mean <- mean(degrees, na.rm = T)
  degree_maximum <- max(degrees, na.rm = T)
  ## Density
  # Density: The proportion of present edges from all possible edges in the network.
  density <- igraph::edge_density(graph)
  ## Diameter
  # Diameter: the longest geodesic distance (length of the shortest path between two nodes) in the network. In igraph, diameter() returns the distance, while get_diameter() returns the nodes along the first found path of that distance.
  diameter <- igraph::diameter(graph, directed = F, unconnected = TRUE, weights = (1 / E(graph)$weight)) ### here use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
  ## Efficiency
  # Network efficiency: average connection efficiency of all pairs of nodes, where connection efficiency is the reciprocal of the shortest path length between the two nodes
  net_dist <- shortest.paths(graph, weights = 1 / E(graph)$weight, mode = "all") ## again use the inverse of the weights, because the algorithm considers weights as distances rather than strengths of connexion
  net_dist[net_dist == 0] <- NA ## remove distances to self
  efficiency <- 1 / net_dist ## transform each distance into an efficiency
  efficiency <- (1 / ((vcount(graph) * (vcount(graph) - 1)))) * (sum(efficiency, na.rm = TRUE))
  ## Modularity
  communities <- cluster_louvain(graph, weights = E(graph)$weight)
  community_membership <- communities$membership
  modularity <- modularity(graph, community_membership, weights = E(graph)$weight)

  # FINAL OUTPUT #DATAFRAME with the network properties per each 3 hours timeslot
  summary_collective <- rbind(summary_collective, data.frame(
    task_assortativity = task_assortativity,
    clustering = clustering,
    degree_mean = degree_mean,
    degree_maximum = degree_maximum,
    density = density,
    diameter = diameter,
    efficiency = efficiency,
    modularity = modularity, stringsAsFactors = F
  ))

  ##### INDIVIDUAL NETWORK PROPERTIES ######################################
  #### inherited from Stroeymeyt et al. 2018
  #### Part 2: individual network properties ####

  #### Part 2: individual network properties ####
  ### prepare table
  individual <- data.frame(
    antID = V(graph)$antID,
    IsQueen = V(graph)$IsQueen,
    Exposed = V(graph)$Exposed,
    AntTask = V(graph)$AntTask,
    AntTask_num = V(graph)$AntTask_num,
    degree = "NULL",
    aggregated_distance_to_queen = "NULL" # ,
    # mean_aggregated_distance_to_treated=NA,
    # same_community_as_queen=NA
  )
  ## degree
  # individual$degree <- degrees
  individual[match(names(degrees), individual$antID), "degree"] <- degrees
  # ##same community as queen
  # queen_comm <- community_membership[which(V(net)$name==queenid)]
  # community_membership <- community_membership==queen_comm
  # individual[match(V(net)$name,individual$tag),"same_community_as_queen"] <- community_membership
  ## path length to queen
  path_length_to_queen <- t(shortest.paths(graph, v = V(graph), to = V(graph)[which(V(graph)$IsQueen == TRUE)], weights = 1 / E(graph)$weight))
  individual[match(colnames(path_length_to_queen), individual$antID), "aggregated_distance_to_queen"] <- as.numeric(path_length_to_queen)
  ######## Mean path length to treated; aggregated_network
  path_length_to_treated <- as.data.frame(as.matrix(shortest.paths(graph, v = V(graph), to = V(graph)[which(V(graph)$Exposed == TRUE)], weights = 1 / E(graph)$weight)))
  path_length_to_treated["mean_distance_to_treated"] <- NA
  path_length_to_treated$mean_distance_to_treated <- as.numeric(rowMeans(path_length_to_treated, na.rm = T))
  individual[match(rownames(path_length_to_treated), individual$antID), "mean_aggregated_distance_to_treated"] <- path_length_to_treated[, "mean_distance_to_treated"]

  # --------
  # Get complete list of Ants
  AntID_list <- NULL
  for (ant in 1:length(e$ants)) {
    AntID_list <- c(AntID_list, e$ants[[ant]]$ID)
  }

  # #add missing ants
  missing_ants <- subset(AntID_list, !(AntID_list %in% individual$antID))
  missing_ants_table <- data.frame()

  for (MISSING in missing_ants) {
    for (id in length(e$ants[[MISSING]]$identifications)) {
      # print ( e$ants[[MISSING]]$identifications[[id]]$tagValue )
      missing_ants_table <- rbind(missing_ants_table, data.frame(antID = MISSING))
    }
  }

  if (nrow(missing_ants_table) > 0) {
    # add empty cols
    missing_ants_table[setdiff(names(individual), names(missing_ants_table))] <- NA
    individual <- rbind(individual, missing_ants_table)
  }
  individual <- individual[order(individual$antID), ]

  # -------

  ### Add data to main data table
  summary_individual <- rbind(summary_individual, individual) ### NOT SURE THAT THIS STEP WORKS!!! AS AFTER THE FIRST RUN IT BECOMES A LIST OBKJECT AND RBIND MAY MESS0UP THINGS


  ####################
  # return a list object.

  summaries_network <- list(summary_collective = summary_collective, summary_individual = summary_individual)
  return(summaries_network)
}
