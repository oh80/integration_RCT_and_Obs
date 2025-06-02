add_meta_info <- function(data, ...){
  data$info = list(...)
  return(data)
}

count_num_file <- function(Date_path, file_name){
  file_list <- list.files(Date_path)
  
  same_name_file <- file_list[base::grepl(pattern = file_name, file_list)]
  
  output <- length(same_name_file) + 1
  return(output)
}

make_save_path <- function(data_type, sample_size){
  # make date folder
  Date <- Sys.time() |> format("%m%d ") |> trimws()
  Date_path <- here::here("01_data", "data", Date)
  if (!dir.exists(Date_path)) {
    dir.create(Date_path, recursive = TRUE)
  }

  file_name <- paste0(data_type, "_n" , sample_size)
  file_num = count_num_file(Date_path, file_name) 
  
  file_path <- paste0(file_name, "_", file_num, ".obj")
  path <- here::here("01_data", "data", Date, file_path)

  return(path)
}


save_data <- function(data, data_type, sample_size){
  path <- make_save_path(data_type, sample_size)
  saveRDS(data, path)
}