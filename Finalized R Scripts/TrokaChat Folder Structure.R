#' TrokaChat Folder Structure
#'
#' This function creates a specific folder structure within a specified directory.
#' The main folder is named "TrokaChat", and within this folder, it creates
#' several subfolders ("csv1", "TrokaChat Output", "TrokaChat Compare",
#' "TrokaChat Charts", "For GSEA", "Tensor", "GSEA", "DEG Output").
#' Additionally, within the "TrokaChat Charts" folder, it creates three
#' sub-subfolders ("save", "save2", "save3").
#'
#' @param main_directory A character string specifying the path to the main directory
#' under which the folder structure will be created.
#'
#' @return This function does not return a value. It  is called for its side effect,
#' which is creating a specific folder structure within the specified directory.
#'
#' @examples
#' \dontrun{
#' TrokaChat_folder_structure("/Users/Troks27/Desktop/8. Leticia All/Python TrokaChat")
#' }
#'
#' @export
TrokaChat_folder_structure <- function(main_directory) {
  # Define the folder names
  folders <- c("csv1", "TrokaChat Output", "TrokaChat Compare", "TrokaChat Charts",
               "For GSEA", "Tensor", "GSEA", "DEG Output")

  subfolders <- c("save", "save2", "save3")

  # Create the main folder
  main_folder <- file.path(main_directory, "TrokaChat")
  if (!dir.exists(main_folder)) {
    dir.create(main_folder)
  }

  # Create the subdirectories in the main folder
  for (folder in folders) {
    subdirectory <- file.path(main_folder, folder)
    if (!dir.exists(subdirectory)) {
      dir.create(subdirectory)
    }
  }

  # Create the sub-subdirectories in the TrokaChat Charts folder
  for (subfolder in subfolders) {
    sub_subdirectory <- file.path(main_folder, "TrokaChat Charts", subfolder)
    if (!dir.exists(sub_subdirectory)) {
      dir.create(sub_subdirectory)
    }
  }
}
