add_meta_info <- function(data, ...){
  data$info = list(...)
  return(data)
}


count_num_file <- function(Date, data_type){
  path <- here::here("01_data", "data")
  file_list <- list.files(path)
  
  same_Date_file <- file_list[base::grepl(pattern = Date, file_list)]
  same_type_file <- same_Date_file[base::grepl(pattern = data_type, same_Date_file)]
  
  output <- length(same_type_file) 
  return(output)
}


make_save_path <- function(data_type){
  Date <- Sys.time() |> format("%m%d ") |> trimws()
  data_number = count_num_file(Date, data_type) + 1

  file_path <- paste0(data_type, "_", Date, "_", data_number, ".obj")
  path <- here::here("01_data", "data", file_path)
  return(path)
}


save_data <- function(data, data_type){
  path <- make_save_path(data_type)
  saveRDS(data, path)
}