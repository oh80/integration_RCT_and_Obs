source(here::here("04_report", "code", "utils.R"))

main <- function(){
  #setting
  analyze_date <- "1010"
  analyze_name <- "proposal_1d_squared_n250_2_1.obj"
  method <- "proposal"
  #plot_params <- c("eta_g", "eta_tau", "eta_b")
  plot_params <- c("mu_b", "l_tau", "l_b")
  
  # MCMC samples
  analyze_path <- here::here("03_analyze", "result", analyze_date, analyze_name)
  analyze_res <- readRDS(analyze_path)
  
  # plot
  trace_plot <- plot_sample_path(analyze_res$samples, plot_params)
  save_plot(trace_plot, paste0(analyze_name,"_trace_plot"))
}


plot_sample_path <- function(samples, params){
  plot_list <- list()
  for(param in params){
    df_for_plot <- data.frame("t"=seq(1, dim(samples[[param]])[1], 1),
                              param = samples[[param]])
    plot <- ggplot2::ggplot(data=df_for_plot,
                            mapping = ggplot2::aes(x=t, y=param)) +
      ggplot2::geom_line()+ 
      ggplot2::theme_gray()  +
      ggplot2::ylim(c(0, 10)) + 
      ggplot2::labs(y=param)
    
   plot_list[[param]] = plot
  }
  joined_plot <- gridExtra::grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]])
  return(joined_plot)
}

main()