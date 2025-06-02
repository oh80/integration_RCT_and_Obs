extract_file_name <- function(path){
  file_name <- sub(".*01_data/data/\\d{4}/", "", path)
  file_name <- sub(".obj", "", file_name)
  return(file_name)
}


count_num_file <- function(Date_path, file_name){
  file_list <- list.files(Date_path)
  
  same_name_file <- file_list[base::grepl(pattern = file_name, file_list)]
  
  output <- length(same_name_file) + 1
  return(output)
}


make_res_path <- function(model, data_path){
  # make date folder
  Date <- Sys.time() |> format("%m%d ") |> trimws()
  Date_path <- here::here("03_analyze", "result", Date)
  
  if (!dir.exists(Date_path)) {
    dir.create(Date_path, recursive = TRUE)
  }
  
  # make save path
  data_name <- extract_file_name(data_path)
  file_name <- paste0(model, "_" , data_name)
  file_num = count_num_file(Date_path, file_name) 
  
  file_path <- paste0(file_name,"_", file_num, ".obj")
  path <- here::here("03_analyze", "result", Date, file_path)
}


save_result <- function(result, model, data_path){
  path <- make_res_path(model, data_path)
  saveRDS(result, path)
}