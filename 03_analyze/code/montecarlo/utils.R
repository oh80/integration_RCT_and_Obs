count_num_file <- function(Date_path, file_name){
  file_list <- list.files(Date_path)
  
  same_name_file <- file_list[base::grepl(pattern = file_name, file_list)]
  
  output <- length(same_name_file) + 1
  return(output)
}


make_res_path <- function(setting_name){
  # make date folder
  Date <- Sys.time() |> format("%m%d ") |> trimws()
  Date_path <- here::here("03_analyze", "result", Date)
  
  if (!dir.exists(Date_path)) {
    dir.create(Date_path, recursive = TRUE)
  }
  
  # make save path
  file_name <- paste0(setting_name, "_montecarlo")
  file_num  <-  count_num_file(Date_path, file_name) 
  
  file_path <- paste0(file_name,"_", file_num, ".obj")
  path <- here::here("03_analyze", "result", Date, file_path)
}


save_result <- function(result, setting_name){
  path <- make_res_path(setting_name)
  saveRDS(result, path)
}