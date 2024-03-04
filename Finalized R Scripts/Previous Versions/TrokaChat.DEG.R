#' TrokaChat.DEG function
#'
#' @description
#' This function conducts a multi-step differential expression gene (DEG) analysis with filtering and cluster analysis.
#'
#' @param object Seurat object containing gene expression data.
#' @param samples Numeric, indicating the number of samples in the study.
#' @param shortidents Vector, sample identifiers, control condition should be listed first.
#' @param filepath String, file path for saving the outputs.
#' @param export_folder_1 String, name of the folder where the CSV outputs should be saved.
#' @param clusters String, name of the column in the Seurat object metadata containing cluster assignments.
#' @param sample_clus String, name of the column in the Seurat object metadata containing the sample cluster assignments.
#' @param cluster_range Numeric vector, range of clusters to consider in the analysis.
#' @param sample_species String, species of the samples ("human" or "mouse").
#' @param cluster_pct_thresh Numeric, threshold for cluster percentage (between 0 and 100).
#' @param directory_path String, path to the directory containing the signaling pathways database.
#'
#' @return The function saves CSV files containing DEG results for each group to the specified folder and does not return any value.
#'
#' @examples
#' TrokaChat.DEG(object = merged_seurat_all_finalsubset_filtered,
#'               samples = 3,
#'               shortidents = c("H","NL","LS"),
#'               filepath = "./TrokaChat/",
#'               export_folder_1 = "csv1",
#'               clusters = "combined_cluster_number",
#'               sample_clus = "sample_clus",
#'               cluster_range = 0:11,
#'               sample_species = "human",
#'               cluster_pct_thresh = 0,
#'               directory_path = "/Users/Troks27/Desktop/8. Leticia All/Python TrokaChat/")
#'
#' @export
#' @importFrom Seurat Idents SetIdent FindMarkers DefaultAssay
#' @importFrom utils write.csv
#' @importFrom base factor paste0 stop rm
#' @importFrom gtools mixedsort
#' @importFrom openxlsx read.xlsx
#' @importFrom stats na.omit
#' @importFrom stringr str_detect str_replace
#' @importFrom purrr map_chr
#' @importFrom grDevices col2rgb

TrokaChat.DEG <- function(object, samples, shortidents, filepath, export_folder_1,
                          clusters, sample_clus, cluster_range, sample_species,
                          cluster_pct_thresh, directory_path) {
  ## Selecting Correct Ligand-Receptor (LR) Database for Species based on sample_species ----
  # Then, it loads the respective database and creates a vector of signaling genes.
  # Input: Seurat object, sample_species ("human" or "mouse")
  # Output: all_signaling_genes_final
  Idents(object) <- factor(x = Idents(object), levels = mixedsort(levels(object)))
  Idents(object)

  # Construct the file path based on sample_species
  if (sample_species == "human") {
    file_name <- "human_ALL_PATHWAYS_UPDATED.xlsx"
  } else if (sample_species == "mouse") {
    file_name <- "All Signaling Pathways.xlsx"
  } else {
    # Handle the case when sample_species does not match either condition
    stop("Invalid sample_species value. Please provide 'human' or 'mouse'.")
  }

  # Construct the full file path
  file_path <- paste0(directory_path, file_name)

  # Read the Excel file
  lr_database <- read.xlsx(file_path)

  # Extract values from the 'ligand' column
  ligand_values <- as.character(lr_database$ligand)

  # Extract values from the 'subunit_' columns
  subunit_columns <- grep("^subunit_", names(lr_database), value = TRUE)
  subunit_values <- unlist(lr_database[subunit_columns])
  subunit_values <- as.character(subunit_values)

  # Combine the values into a single vector
  signaling_genes <- c(ligand_values, subunit_values)
  all_signaling_genes_final <- signaling_genes


  ## Generating Counts Table based on clusters and shortidents. ----
  # Creates a table that shows the number of cells per cluster for each condition.
  # Input: Seurat object, clusters, shortidents
  # Output: n_cells_bycondition_percluster
  object <- SetIdent(object, value = clusters)
  Idents(object) <- factor(x = Idents(object), levels = mixedsort(levels(object)))
  #Idents(object)
  n_cells_bycondition_percluster<-table(object@active.ident, object@meta.data$sample)
  n_cells_bycondition_percluster <- n_cells_bycondition_percluster[ , shortidents]
  #n_cells_bycondition_percluster


  ## Create List of Cluster Combinations for DEG analysis ----
  # Creates a list of all possible cluster combinations for subsequent DEG analysis
  # Input: shortidents, cluster_range, n_cells_bycondition_percluster
  # Output: total_list
  group <- paste0(shortidents, "_vs_", shortidents[1])

  total_list <- vector("list", length(shortidents))
  list_of_lists <- vector("list", length(shortidents))


  for (i in seq_along(shortidents)) {
    my_list <- vector("list", length(shortidents))

    cutoff <- sum(n_cells_bycondition_percluster[, i]) * (cluster_pct_thresh / 100)

    new_element <- names(which(n_cells_bycondition_percluster[, i] <= 3 |
                                 n_cells_bycondition_percluster[, i] <= floor(cutoff)))
    my_list[[i]] <- if(length(new_element) == 0) NA else new_element

    my_list_i <- paste0(shortidents[i], "_", my_list[[i]])

    my_list_2_i = {}
    # Create the list using nested lapply
    my_list_2_i <- lapply(1:length(shortidents), function(i) {
      ident_values <- lapply(cluster_range, function(a) {
        ident_name <- shortidents[i]
        paste0(ident_name, "_", a)
      })
      # Flatten the list to a vector and sort it
      mixedsort(unlist(ident_values))
    })

    my_list_3_i <- setdiff(my_list_2_i[[1]], my_list_i)

    clus_nums <- as.numeric(gsub(".*?([0-9]+).*", "\\1", my_list_3_i))
    my_list_4_i <- as.list(clus_nums)

    list_i <- list(unlist(my_list_4_i))

    list_i_final <- list_i
    list_of_lists[[i]] <- list_i[[1]]

    for (a in seq_along(my_list_4_i)) {
      rem <- my_list_4_i[[a]]
      element <- list_i[[1]][!list_i[[1]] %in% rem]
      element <- purrr::map_chr(element, ~ paste0(shortidents[1], "_", .))
      list_i_final[[a]] <- element
    }

    labels <- purrr::map_chr(my_list_4_i, ~ paste0(shortidents[i], "_", .))
    names(list_i_final) <- labels

    total_list[[i]] <- list_i_final
  }
  print(total_list)
  saveRDS(total_list, file = "./TrokaChat/csv1/total_list.rds")



  ##Initial DEG (Differential Expression Gene) Analysis ----
  # Executes initial DEG analysis by running FindMarkers() for each cluster combination in total_list
  # Input: Seurat object, total_list, all_signaling_genes_final
  # Output: Objects named *_TrokaChat1, *_TrokaChat2, *_TrokaChat3, etc. in the global environment
  rm(list=ls(pattern="_TrokaChat1"))
  rm(list=ls(pattern="_TrokaChat2"))
  rm(list=ls(pattern="_TrokaChat3"))

  object <- SetIdent(object, value = sample_clus)
  DefaultAssay(object) <- "RNA"
  features = {}
  for (i in 1:length(total_list)){
    for(a in 1:length(total_list[[i]]))
      assign(print(paste0(group[i],paste0(gsub("\\D+", "", names(total_list[[i]][a]))),"_TrokaChat",i)), FindMarkers(object, logfc.threshold = 0, only.pos = FALSE, ident.1 = print(paste0(names(total_list[[i]][a]))), ident.2 = print(total_list[[i]][[a]]), features=intersect(rownames(object), all_signaling_genes_final)), envir = environment())
  }




  ## For each cluster, only those DEGs that have a p-value less than 0.05 and are present in at least one of the conditions tested are considered significant. These genes are then used in subsequent DEG analysis. ----
  # Input: DEG results from the previous step
  # Output: TC_overall_list
  # Create the datalist
  datalist <- list()

  for (i in seq_along(shortidents)) {
    pattern_matches <- grep(paste0("*._TrokaChat",i,"$"), ls(environment()), value = TRUE)
    pattern <- mixedsort(pattern_matches)
    pattern_list <- do.call("list", mget(pattern))

    # str(pattern_list)  # Uncomment if you want to print the structure of pattern_list

    pattern <- gsub(paste0("_TrokaChat",i), '', pattern)
    pattern <- gsub('\\D', '', pattern)

    for (a in seq_along(pattern_list)) {
      if (nrow(pattern_list[[a]]) > 0) {
        pattern_list[[a]]$cluster <- pattern[a]
      }
    }

    datalist[[i]] <- do.call(rbind, pattern_list)
  }

  TC_all <- do.call(rbind, datalist)

  print(head(TC_all))
  # Clean up the gene names and filter
  TC_all$gene <- gsub("^.*\\.","", rownames(TC_all))
  TC_all <- TC_all[TC_all$p_val_adj < 0.05,]
  TC_all <- TC_all[TC_all$gene %in% all_signaling_genes_final, ]
  TC_all <- TC_all[order(TC_all$cluster),]

  # Split, remove duplicates, and sort
  TC_all_split <- split(TC_all, f = TC_all$cluster)
  TC_all_split <- lapply(TC_all_split, function(x) x[!duplicated(x[c("gene")]), ])
  TC_all_split <- TC_all_split[mixedsort(names(TC_all_split))]

  # Print the structure
  str(TC_all_split)

  # Create the TC_overall_list
  TC_overall_list<- lapply(1:length(total_list), function(i) {
    temp <- TC_all_split
    list <- setdiff(names(temp), list_of_lists[[i]])
    temp[names(temp) %in% list == FALSE]
  })
  ## Conduct DEG With filtered list of genes ----
  # Run DEG analysis again but now using only significant genes from TC_overall_list.
  # Input: Seurat object, total_list, TC_overall_list
  # Output: Objects named *_TrokaChat_a, *_TrokaChat_b, *_TrokaChat_c, etc. in the global environment
  rm(list=ls(pattern="_TrokaChat_"))

  for (i in 1:length(total_list)){
    for (a in 1:length(TC_overall_list[[i]])){
      assign(print(paste0(group[i],paste0(gsub("\\D+", "", names(total_list[[i]][a]))),"_TrokaChat_",letters[i])), FindMarkers(object, logfc.threshold = 0,features = TC_overall_list[[i]][[a]]$gene, only.pos = FALSE, ident.1 = print(paste0(names(total_list[[i]][a]))), ident.2 = print(total_list[[i]][[a]])), envir = environment())
    }
  }

  saveRDS(TC_overall_list, file = "./TrokaChat/csv1/TC_overall_list.rds")
  ## Assembles DEG data and pulls genes for final DEG anaylsis to recover pct.2 ----
  # Aggregates the DEG analysis results and prepares for the final DEG analysis.
  # Input: Results from the second DEG analysis
  # Output: dataframe_names_final
  rm(list=ls(pattern="dataframe_names_"))

  # Preallocate a list to hold your dataframes
  dataframe_names_list <- vector("list", length(total_list))
  Pattern2_list <- vector("list", length(total_list))

  # First loop
  for (i in seq_along(total_list)) {
    dataframe_names_list[[i]] <- mget(mixedsort(ls(pattern = paste0("_TrokaChat_", letters[i]))))
  }

  # Second loop
  for (i in seq_along(total_list)) {
    Pattern2_list[[i]] <- gsub(paste0("_TrokaChat_", letters[i]), '', names(dataframe_names_list[[i]]))
    Pattern2_list[[i]] <- gsub('\\D', '', Pattern2_list[[i]])
    perm2 <- Pattern2_list[[i]]
    var <- dataframe_names_list[[i]]

    for (a in seq_along(perm2)) {
      var[[a]]$cluster <- perm2[[a]]
      var[[a]]$gene <- row.names(var[[a]])
    }
    dataframe_names_list[[i]] <- var
  }



  ## Run Final DEG ----
  # Perform the final DEG analysis using the filtered list of genes.
  # Input: Seurat object, total_list, dataframe_names_list
  # Output: dataframe_names_new_list
  rm(list=ls(pattern="TrokaChatnew"))
  rm(list=ls(pattern="_new"))

  # Create a list to hold new dataframes
  dataframe_names_new_list <- vector("list", length(total_list))

  # Outer loop
  for (i in seq_along(total_list)){

    # Initialize a list to temporarily hold the FindMarkers results for each 'a'
    FindMarkers_results <- vector("list", length(dataframe_names_list[[i]]))

    # Inner loop
    for (a in seq_along(dataframe_names_list[[i]])){
      holder <- dataframe_names_list[[i]]

      # Create the new variable name
      new_var_name <- paste0(group[i], paste0(gsub("\\D+", "", names(total_list[[i]][a]))), "_TrokaChatnew_", letters[i])
      print(new_var_name)
      print(names(total_list[[i]][a]))
      print(names(total_list[[i]][-a]))

      # Call FindMarkers and store the result in the list
      FindMarkers_results[[a]] <- FindMarkers(
        object,
        features = holder[[a]]$gene,
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = FALSE,
        ident.1 = names(total_list[[i]][a]),
        ident.2 = names(total_list[[i]][-a])
      )
    }
    # Collect all FindMarkers_results for this 'i' and store in dataframe_names_new_list
    dataframe_names_new_list[[i]] <- FindMarkers_results
  }


  ## Add the pct.2 values from this analysis back to DEG output ----
  # Replaces the pct.2 column in dataframe_names_list with the one from dataframe_names_new_list
  # Input: dataframe_names_new_list, dataframe_names_list
  # Output: dataframe_names_final
  # Loop over each 'a' in total_list
  for (a in seq_along(total_list)) {

    # Loop over each 'i' in dataframe_names_new_list[[a]]
    for (i in seq_along(dataframe_names_new_list[[a]])) {
      temp <- dataframe_names_new_list[[a]]
      temp2 <- dataframe_names_list[[a]]

      # Order rows by row.names
      temp[[i]] <- temp[[i]][order(row.names(temp[[i]])), ]
      temp2[[i]] <- temp2[[i]][order(row.names(temp2[[i]])), ]

      # Replace pct.2 column in temp2 with the one from temp
      temp2[[i]]$pct.2 <- temp[[i]]$pct.2

      # Update dataframe_names_list
      dataframe_names_list[[a]] <- temp2
    }
  }

  # Initialize a list to hold the final dataframes
  dataframe_names_final <- vector("list", length(total_list))

  # Loop over each 'i' in total_list
  for (i in seq_along(total_list)) {
    # Concatenate dataframes together and store the result in dataframe_names_final
    dataframe_names_final[[i]] <- do.call(rbind, dataframe_names_list[[i]])
  }


  ## Export the csv outputs ----
  # Writes the final DEG results to CSV files, one for each group
  # Input: dataframe_names_final, filepath, export_folder_1, group
  # Output: CSV files
  for(i in seq_along(dataframe_names_final)){
    write.csv(dataframe_names_final[[i]],file=paste0(filepath,export_folder_1,"/",paste0(group[i]),".csv"))
  }

}
