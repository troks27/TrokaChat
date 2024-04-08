#' Perform TrokaChat Null Distribution Computation
#'
#' This function performs the TrokaChat null distribution computation using the provided arguments.
#'
#' @param object The Seurat object to be used for the analysis.
#' @param shortidents A character vector representing short identifiers for conditions. The control condition should be listed first.
#' @param filepath The directory path where the output files should be saved.
#' @param export_folder The name of the export folder within the filepath.
#' @param clusters The name of the column in the Seurat object representing the clusters.
#' @param n.perms The number of permutations to perform.
#'
#' @return CSV files containing the null distribution DEGs for each condition. These files are stored in the specified export folder.
#'
#' @examples
#' TrokaChat.DEG.nulldist(object = merged_seurat_all_finalsubset_filtered,
#'                        shortidents = c("H","NL","LS"),
#'                        filepath = "./TrokaChat/",
#'                        export_folder = "csv1",
#'                        clusters = "combined_cluster_number",
#'                        n.perms = 50)
#'
#' @importFrom Seurat SetIdent FetchData RenameCells
#' @importFrom dplyr count
#' @importFrom tidyr spread
#' @importFrom data.table rbindlist
#' @export
TrokaChat.DEG.nulldist <- function(object, shortidents, filepath, export_folder, clusters, n.perms) {

  export_folder_1 = export_folder
  group <- paste0(shortidents, "_vs_", shortidents[1])

  ## Create counts table ----
  object <- SetIdent(object, value = paste0(clusters))

  n_cells_bycondition_percluster <- FetchData(object, vars = c("ident", "sample")) %>%
    dplyr::count(ident, sample) %>%
    tidyr::spread(ident, n)

  filename <- file.path(filepath, export_folder_1, "counts.csv")

  write.csv(n_cells_bycondition_percluster, file = filename, row.names = FALSE)

  total_list <- readRDS(file = "./TrokaChat/csv1/total_list.rds")
  TC_overall_list <- readRDS(file = "./TrokaChat/csv1/TC_overall_list.rds")
  ## creating duplicate object----
  holder2 <- list()
  object_new<-object
  object_new@meta.data$TrokaChat_clusters_new <- object_new@meta.data[[clusters]]
  ## rename_and_merge_object function ----
  rename_and_merge_object <- function(object_new, shortidents){
    object_sub <- subset(x = object_new, subset = sample == shortidents[1])
    object_sub2 <- subset(x = object_new, subset = sample == shortidents[1])
    object_sub <- RenameCells(object_sub, c('cell1'))
    object_sub@meta.data$sample <- paste0(shortidents[1],'_NEW')
    object_final <- merge(object_sub, object_sub2)
    return(object_final)
  }
  ## Find Markers and Create DF function ----
  findMarkersAndCreateDF <- function(object, total_list, group, c, TC_overall_list) {
    # Initialize a list to store the data frames directly
    dataframe_names_fornull <- list()

    # Initialize variables to store unique ident.1 and ident.2 values for confirmation
    unique_ident1 <- unique(sapply(total_list[[c]], names))
    unique_ident2_lengths <- integer(length(unique_ident1))

    for (i in 1:length(TC_overall_list)){
      # Determine ident.1 and ident.2 for the current iteration
      ident1 <- names(total_list[[c]][i])
      ident2 <- total_list[[c]][[i]]

      # Perform marker finding and store directly in the list
      dataframe_names_fornull[[i]] <- FindMarkers(object,
                                                  features = TC_overall_list[[i]]$gene,
                                                  logfc.threshold = 0,
                                                  min.pct = 0,
                                                  only.pos = FALSE,
                                                  ident.1 = ident1,
                                                  ident.2 = ident2
      )
      # Assign gene names as a new column in each dataframe
      dataframe_names_fornull[[i]]$gene <- row.names(dataframe_names_fornull[[i]])

      # Extract the cluster number from ident1 after the "_"
      cluster_num <- sub(".*_", "", ident1) # This will capture the part after the last underscore in ident1
      dataframe_names_fornull[[i]]$cluster <- as.numeric(cluster_num) # Ensure it's numeric

      # Update ident.2 lengths for confirmation
      unique_ident2_lengths[which(unique_ident1 == ident1)] <- length(ident2)
    }

    # Print a summary confirmation message
    confirmation_message <- paste("All permutations processed. Confirmed ident.1 as condition 1_cluster for", length(unique_ident1), "permutations. Each ident.2 correctly excluded the ident.1 cluster and included all other clusters from condition 2.")
    print(confirmation_message)

    # Return the list of data frames directly
    return(dataframe_names_fornull)
  }

  ## Process Dataframes function ----
  processDataframes <- function(df_list) {
    for( i in seq_along(df_list)){
      df_list[[i]]$gene <-  row.names(df_list[[i]])
    }
    list2env(df_list,envir=.GlobalEnv)

    Pattern3 <- df_list
    nombres <- names(Pattern3)
    Pattern3 <- gsub(paste0("_TrokaChat_nulldist"), '', nombres)
    Pattern3 <- gsub('\\D', '', Pattern3)
    perm2 <- Pattern3
    var <- df_list

    for (b in seq_along(perm2)){
      var[[b]]$cluster <- perm2[[b]]
    }
    df_list <- var
    return(df_list)
  }
  ## Find Markers second function ----
  FindMarkersagain <- function(total_list, object, c, TC_overall_list, dataframe_names_fornull) {
    dataframe_names_fornull_new <- list()

    for (i in seq_along(TC_overall_list)) {

      ident1 <- names(total_list[[c]][i])
      ident2 <- names(total_list[[c]][-i])

      # Perform marker finding using the correct features source
      result <- FindMarkers(
        object,
        features = dataframe_names_fornull[[i]]$gene, # Corrected features source
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = FALSE,
        ident.1 = ident1,
        ident.2 = ident2
      )

      # Assign gene names as a new column in the result
      result$gene <- row.names(result)

      # Extract the cluster number from ident1 after the "_"
      cluster_num <- sub(".*_", "", ident1) # Extract cluster number from ident1
      result$cluster <- as.numeric(cluster_num) # Ensure it's numeric

      dataframe_names_fornull_new[[i]] <- result
    }
    return(dataframe_names_fornull_new)
  }






  ## Process New Dataframes function ----
  process_data <- function(a, holder2, dataframe_names_fornull, dataframe_names_fornull_new) {
    for (i in seq_along(dataframe_names_fornull_new)){

      # Assuming 'row.names' can be converted to a column for joining purposes
      df_fornull <- dataframe_names_fornull[[i]] %>%
        tibble::rownames_to_column("row_name")
      df_fornull_new <- dataframe_names_fornull_new[[i]] %>%
        tibble::rownames_to_column("row_name")

      df_fornull <- df_fornull[order(df_fornull$gene), ]
      df_fornull_new <- df_fornull_new[order(df_fornull_new$gene), ]



      #print("Head of df_fornull")
      #print(head(df_fornull))
      #print("Head of df_fornull_new")
      #print(head(df_fornull_new))

      # Perform the left join and update pct.2 directly
      df_fornull_updated <- df_fornull %>%
        left_join(df_fornull_new %>% select(row_name, `pct.2`), by = "row_name", suffix = c("", "_new")) %>%
        mutate(`pct.2` = coalesce(`pct.2_new`, `pct.2`)) %>%
        select(-`pct.2_new`) %>%
        tibble::column_to_rownames("row_name")


      #print("Head of df_fornull_updated")
      #print(head(df_fornull_updated))


      dataframe_names_fornull[[i]] <- df_fornull_updated

      #print("Head of dataframe_names_fornull[[i]]")
      #print(head(dataframe_names_fornull[[i]]))
    }

    dataframe_names_fornull_final <- do.call(rbind, dataframe_names_fornull)
    holder2[[a]] <- dataframe_names_fornull_final
    print(a)
    print("This has corrected pct.2")


    # Iterate over each data frame stored in holder2
    holder2 <- lapply(holder2, function(df) {
      # Ensure that both 'cluster' and 'gene' columns exist
      if("cluster" %in% names(df) && "gene" %in% names(df)) {
        # Sort the data frame by 'cluster' and then by 'gene' alphabetically
        df_sorted <- df %>%
          arrange(cluster, gene)
        return(df_sorted)
      } else {
        # Return the original data frame if 'cluster' or 'gene' column is missing
        return(df)
      }
    })

    #str(holder2)

    rm(list=ls(pattern="_TrokaChat_newnull"))
    rm(list=ls(pattern="_TrokaChat_nulldist"))
    rm(list=ls(pattern="dataframe_names_fornull"))
    rm(list=ls(pattern="dataframe_names_fornull_new"))
    rm(list=ls(pattern="dataframe_names_fornull_final"))

    return(holder2)
  }
  ## Process holder Function ----
  process_holder <- function(export_folder, group_element) {
    holder3 <- NULL
    holder3 <- rbindlist(holder2)[,lapply(.SD,mean), list(gene, cluster)]
    holder3 <- as.data.frame(holder3)
    col_order <- c("p_val", "avg_log2FC", "pct.1","pct.2", "p_val_adj", "gene", "cluster")
    holder3 <- holder3[, col_order]

    # Create a new folder for all of the output files from one comparison, and set the filepath below.
    write.csv(holder3, file = paste0(filepath, export_folder, "/nulldistavgs_", group_element, ".csv"), row.names = FALSE)
  }
  ## Main Script ----
  for (c in 1:length(total_list)){
    if (c==1){
      object_temp <- object_new
      object_final <- rename_and_merge_object(object_temp, shortidents)


      # Step 1: Check the Active Assay
      activeAssay <- DefaultAssay(object_final)

      # Step 2: Conditional Execution
      if(activeAssay == "SCT") {
        # Step 3: Apply Function
        object_final <- PrepSCTFindMarkers(object = object_final)
        # Optionally, print a message indicating the operation was performed
        message("PrepSCTFindMarkers has been applied to the object.")
      } else {
        # Optional: Message to indicate the active assay is not SCT, hence no action taken
        message("Active assay is not SCT; no action was taken.")
      }


      for (a in 1:n.perms){

        #New Print Statment
        cat(sprintf("\nStarting permutation iteration: %d\n", a))

        # Identifying the subset directly in object_final@meta.data for permutation
        sample_indices <- which(object_final@meta.data$sample == paste0(shortidents[1],"_NEW") )

        cat("Before permutation:\n")
        print(head(object_final@meta.data$TrokaChat_clusters_new[sample_indices]))

        #Permute the labels
        object_final@meta.data$TrokaChat_clusters_new[sample_indices] <- sample(object_final@meta.data$TrokaChat_clusters_new[sample_indices])

        cat("After permutation:\n")
        print(head(object_final@meta.data$TrokaChat_clusters_new[sample_indices]))

        #View(object_final@meta.data)

        object_final@meta.data[['sample_clus_TrokaChat']] <- ifelse(object_final@meta.data$sample == shortidents[1], paste0(object_final@meta.data$sample,"_", object_final@meta.data[[clusters]]), paste0(object_final@meta.data$sample,"_", object_final@meta.data$TrokaChat_clusters_new))
        object_final <- SetIdent(object_final, value = 'sample_clus_TrokaChat')

        total_list_special <- total_list
        labels <- sub("_", "_NEW_", names(total_list_special[[1]]))
        names(total_list_special[[1]]) <- labels

        #Conserved FindMarkers Function
        print("findMarkersAndCreateDF")
        dataframe_names_fornull = findMarkersAndCreateDF(object_final, total_list_special, group, c, TC_overall_list=TC_overall_list)

        #Conserved second find markers Function
        print("FindMarkers Again")
        dataframe_names_fornull_new <- FindMarkersagain(total_list = total_list, object = object_final,c, TC_overall_list=TC_overall_list, dataframe_names_fornull=dataframe_names_fornull)

        print("process_data")
        #Conserved process dataframes again Function
        holder2 <- process_data(a, holder2, dataframe_names_fornull, dataframe_names_fornull_new)

      }
      str(holder2)
      process_holder(export_folder_1, group[c])
    }

    else{
      rm("dataframe_names_fornull")
      rm("holder2")
      rm("dataframe_names_fornull_new")
      holder2 <- list()
      for (a in 1:n.perms){

        cat(sprintf("\nStarting permutation iteration: %d\n", a))

        sample_indices <- which(object_new@meta.data$sample == paste0(shortidents[c]) )

        cat("Before permutation:\n")
        print(head(object_new@meta.data$TrokaChat_clusters_new[sample_indices]))

        #Permute the labels
        object_new@meta.data <- transform(object_new@meta.data, TrokaChat_clusters_new = sample(TrokaChat_clusters_new))

        cat("After permutation:\n")
        print(head(object_new@meta.data$TrokaChat_clusters_new[sample_indices]))

        object_new@meta.data[['sample_clus_TrokaChat']] <- ifelse(object_new@meta.data$sample == shortidents[c], paste0(object_new@meta.data$sample,"_", object_new@meta.data$TrokaChat_clusters_new), paste0(object_new@meta.data$sample,"_", object_new@meta.data[[clusters]]))
        object_new <- SetIdent(object_new, value = 'sample_clus_TrokaChat')

        print("findMarkersAndCreateDF")
        #Function to replace conserved FindMarkers function
        dataframe_names_fornull = findMarkersAndCreateDF(object_new, total_list, group, c, TC_overall_list = TC_overall_list)

        #Conserved second find markers Function
        print("FindMarkers Again")
        dataframe_names_fornull_new <- FindMarkersagain(total_list = total_list, object = object_new,c, TC_overall_list=TC_overall_list, dataframe_names_fornull=dataframe_names_fornull)


        print("process_data")
        #Conserved process dataframes again Function
        holder2 <- process_data(a, holder2, dataframe_names_fornull, dataframe_names_fornull_new)

      }

      str(holder2)
      process_holder(export_folder_1, group[c])

    }
  }
}





