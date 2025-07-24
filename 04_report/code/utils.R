make_res_path <- function(analyze_name, type){
  # make date folder
  Date <- Sys.time() |> format("%m%d ") |> trimws()
  Date_path <- here::here("04_report", "result", Date)
  
  if (!dir.exists(Date_path)) {
    dir.create(Date_path, recursive = TRUE)
  }
  
  # make save path
  
  if(type=="plot"){
    file_name <- paste0(analyze_name,"_plot", ".pdf")
  }else if(type=="metrics"){
    file_name <- paste0(analyze_name, ".obj")
  }
  path <- here::here("04_report", "result", Date, file_name)
  return(path)
}


save_plot <- function(plot, analyze_name){
  path <- make_res_path(analyze_name, type="plot")
  ggplot2::ggsave(path, plot)
}