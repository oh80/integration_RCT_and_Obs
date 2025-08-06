count_num_file <- function(Date_path, file_name){
  file_list <- list.files(Date_path)
  
  same_name_file <- file_list[base::grepl(pattern = file_name, file_list)]
  
  output <- length(same_name_file) + 1
  return(output)
}

make_save_path <- function(name){
  # make date folder
  Date <- Sys.time() |> format("%m%d ") |> trimws()
  Date_path <- here::here("02_build", "data", Date)
  if (!dir.exists(Date_path)) {
    dir.create(Date_path, recursive = TRUE)
  }

  file_num = count_num_file(Date_path, name) 
  
  file_path <- paste0(name, "_", file_num, ".obj")
  path <- here::here("02_build", "data", Date, file_path)
  
  return(path)
}


save <- function(data, name){
  path <- make_save_path(name)
  saveRDS(data, path)
}