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
  ## rename_and_merge_object function ----
  rename_and_merge_object <- function(object_new, shortidents){
    object_sub <- subset(x = object_new, subset = sample == shortidents[1])
    object_sub2 <- subset(x = object_new, subset = sample == shortidents[1])
    object_sub <- RenameCells(object_sub, c('cell1'))
    object_sub@meta.data$sample <- paste0(shortidents[1],'_NEW')
    object_final <- merge(object_sub, object_sub2)
    return(object_final)
  }

  ## creating duplicate object----
  object_new<-object
  holder2 <- list()
  object_new@meta.data$TrokaChat_clusters_new <- object_new@meta.data[[clusters]]


  ## Find Markers and Create DF function ----
  findMarkersAndCreateDF <- function(object, total_list, group, c, TC_overall_list) {
    rm(list=ls(pattern="_TrokaChat_nulldist"))

    for (i in 1:length(TC_overall_list[[c]])){
      assign(print(paste0(group[c],paste0(gsub("\\D+", "", names(total_list[[c]][i]))),"_TrokaChat_nulldist")),
             FindMarkers(object,
                         features = TC_overall_list[[c]][[i]]$gene,
                         logfc.threshold = 0,
                         min.pct = 0,
                         only.pos = FALSE,
                         ident.1 = print(paste0(names(total_list[[c]][i]))),
                         ident.2 = print(total_list[[c]][[i]])
             )
      )
    }
    dataframe_names_fornull <- mget(mixedsort(ls(pattern = "*._TrokaChat_nulldist")))
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
  ## Find Markers scond functin ----
  FindMarkersagain <- function(total_list, object, pattern = "*._TrokaChat_newnull") {
    for (i in seq_along(TC_overall_list[[c]])) {
      assign(
        print(
          paste0(
            group[c],
            paste0(gsub("\\D+", "", names(total_list[[c]][i])), "_TrokaChat_newnull")
          )
        ),
        FindMarkers(
          object,
          features = dataframe_names_fornull[[i]]$gene,
          logfc.threshold = 0,
          min.pct = 0,
          only.pos = FALSE,
          print(paste0(names(total_list[[c]][i]))),
          ident.2 = print(names(total_list[[c]][-i]))
        )
      )
    }
    dataframe_names_fornull_new <- mget(mixedsort(ls(pattern = pattern)))
    return(dataframe_names_fornull_new)
  }
  ## Process New Dataframes function ----
  process_data <- function(a, holder2, dataframe_names_fornull, dataframe_names_fornull_new) {
    for (i in seq_along(dataframe_names_fornull_new)){
      dataframe_names_fornull_new[[i]] <- dataframe_names_fornull_new[[i]][order(row.names(dataframe_names_fornull_new[[i]])), ]
      dataframe_names_fornull[[i]] <- dataframe_names_fornull[[i]][order(row.names(dataframe_names_fornull[[i]])), ]
      dataframe_names_fornull[[i]]['pct.2'] <- dataframe_names_fornull_new[[i]]['pct.2']
    }

    dataframe_names_fornull_final <- do.call(rbind, dataframe_names_fornull)
    holder2[[a]] <- dataframe_names_fornull_final
    print(a)
    str(holder2)

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
      object_final <- rename_and_merge_object(object_new, shortidents)


      for (a in 1:n.perms){
        #Permute the labels
        subset <-subset(object_final@meta.data, sample ==shortidents[1])
        subset <- transform(subset, TrokaChat_clusters_new = sample(TrokaChat_clusters_new))
        inds <- match(object_final@meta.data$Row.names, subset$Row.names)
        object_final@meta.data$TrokaChat_clusters_new[!is.na(inds)] <- subset$TrokaChat_clusters_new[na.omit(inds)]


        object_final@meta.data[['sample_clus_TrokaChat']] <- ifelse(object_final@meta.data$sample == shortidents[1], paste0(object_final@meta.data$sample,"_", object_final@meta.data[[clusters]]), paste0(object_final@meta.data$sample,"_", object_final@meta.data$TrokaChat_clusters_new))
        object_final <- SetIdent(object_final, value = 'sample_clus_TrokaChat')


        total_list_special <- total_list
        labels <- sub("_", "_NEW_", names(total_list_special[[1]]))
        names(total_list_special[[1]]) <- labels

        #Conserved FindMarkers Function
        dataframe_names_fornull = findMarkersAndCreateDF(object_final, total_list_special, group, c, TC_overall_list)

        #Conserved Process Dataframes Function
        dataframe_names_fornull <- processDataframes(dataframe_names_fornull)

        #Conserved second find markers Function
        dataframe_names_fornull_new <- FindMarkersagain(total_list = total_list, object = object_final)

        #Conserved process dataframes again Function
        holder2 <- process_data(a, holder2, dataframe_names_fornull, dataframe_names_fornull_new)

      }
      # Conserved process holder Function
      process_holder(export_folder_1, group[c])
    }

    else{
      for (a in 1:n.perms){
        #Permute the labels
        object_new@meta.data <- transform(object_new@meta.data, TrokaChat_clusters_new = sample(TrokaChat_clusters_new))

        object_new@meta.data[['sample_clus_TrokaChat']] <- ifelse(object_new@meta.data$sample == shortidents[c], paste0(object_new@meta.data$sample,"_", object_new@meta.data$TrokaChat_clusters_new), paste0(object_new@meta.data$sample,"_", object_new@meta.data[[clusters]]))
        object_new <- SetIdent(object_new, value = 'sample_clus_TrokaChat')

        #Function to replace conserved FindMarkers function
        dataframe_names_fornull = findMarkersAndCreateDF(object_new, total_list, group, c, TC_overall_list)

        #Conserved Process Dataframes Function
        dataframe_names_fornull <- processDataframes(dataframe_names_fornull)

        #Conserved second find markers Function
        dataframe_names_fornull_new <- FindMarkersagain(total_list = total_list, object = object_new)

        #Conserved process dataframes again Function
        holder2 <- process_data(a, holder2, dataframe_names_fornull, dataframe_names_fornull_new)

      }
      # Conserved process holder Function
      process_holder(export_folder_1, group[c])
    }
  }
}




