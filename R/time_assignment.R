
#' Evaluate FACS Models for Various Fluorophores
#'
#' This function evaluates models for different fluorophores including BV711, PE, BB515, and BUV737
#' using Generalized Linear Models (GLM).
#' At the end, it returns the filtered dataframe and plots for the GLM fits.
#'
#' @param well_fcs_mc Data frame containing the FACS data with columns for different fluorophores.
#'
#' @return Visualizations for assessing the normality of the fitted model
#' @export
#' @examples
#' # Assuming `my_data` is your FACS data
#' FACS_model_eval(my_data)
#' @section Timebin Assignment
#' @import stats
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import MASS
#' @import Matrix
#' @import dplyr


FACS_model_eval <- function(well_fcs_mc, fluorophores = c('log_BV711.A','log_BUV737.A','log_BB515.A','log_PE.A')) {
  
  #black_list_cells <- which(rowSums(well_fcs_mc[, fluorophores] < 1) > 0)
  black_list_filtered_well_fcs <- well_fcs_mc#[-black_list_cells, ]
  
  fit_list <- list()
  plot_list <- list()
  
  for(fluoro in fluorophores){
    frm <- as.formula(paste0("`",fluoro,"` ~ (`FSC.A` + `FSC.W` + `FSC.H` + `SSC.A` + `SSC.W` + `SSC.H` + `log_APC.A`)^2"))
    glm_fit <- glm(frm,
                   family = gaussian(link = "identity"),
                   data = black_list_filtered_well_fcs)
    
    fit_list[[fluoro]] <- MASS::fitdistr((predict(glm_fit, type = "response") - black_list_filtered_well_fcs[, fluoro]), "normal")
    plot_list[[fluoro]] <- plot(glm_fit, main = paste("GLM fit for", fluoro))
    
  }
    #plot_list <- list()
    #for (i in 1:4) {
    #    #plot(fit, which=i)
    #    plot(glm_fit, main = paste("GLM fit for", fluoro), which = i)
    #    #plot_list[[i]] <- recordPlot()
    #    plot_list[[fluoro]] <- recordPlot()
    #}
    #print(fluoro)
      
    # Set up a 2x2 panel for plots
    #par(mfrow = c(2, 2))
    # Display the plots in the panel
    #options(repr.plot.width=12, repr.plot.height=9)
    #do.call("lapply", list(plot_list, replayPlot))
    # Convert the base plots to grid plots
    #grid_list <- lapply(plot_list, function(p) {
    #  grid.newpage()
    #  grid.echo(function() replayPlot(p))
    #  return(grid.grab())
    #})
    # Display the grid plots in a 2x2 arrangement
    #grid.arrange(grobs = grid_list, ncol = 2)
    #}
    
  #options(repr.plot.width=7, repr.plot.height=7)
  #do.call(grid.arrange, c(plot_list, ncol=2))
  return(plot_list)
}


#' Plot Histograms for FACS Data
#'
#' This function plots histograms for various fluorophores including BV711, PE, BB515, and BUV737.
#' It uses ggplot2 for the plotting.
#'
#' @param well_fcs_mc Data frame containing the FACS data with columns for different fluorophores.
#'
#' @export
#' @examples
#' # Assuming `my_data` is your FACS data
#' plot_FACS_histogram(my_data)
#' @section Timebin Assignment
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import MASS
#' @import Matrix
#' @import dplyr

plot_FACS_histogram <- function(well_fcs_mc, fluorophores = c('log_BV711.A','log_BUV737.A','log_BB515.A','log_PE.A')){
  
  plot_list <- list()
  
  for(fluoro in fluorophores){
    #gg <- ggplot(well_fcs_mc, aes_string(x = paste0('`', fluoro, '`'), fill = 'Stain', colour = 'Stain')) +
    gg <- ggplot(well_fcs_mc, aes_string(x = fluoro, fill = 'Stain', colour = 'Stain')) +
      geom_histogram(aes(y=..density..), breaks=seq(0, 5, 0.05), alpha=0.6, position="identity", lwd=0.2) +zmanr_theme+
      ggtitle(paste("Normalized", fluoro))
    
    plot_list[[fluoro]] <- gg
  }
  
  options(repr.plot.width=12, repr.plot.height=7)
  do.call(ggarrange, c(plot_list, ncol=2, nrow=2))
}

#' FACS_model: Classify cells based on flow cytometry data
#'
#' @param well_fcs_mc Data frame containing flow cytometry data
#' @param fluorophores Vector of fluorophore channels to be considered. Defaults to c('log_BV711-A','log_BUV737-A','log_BB515-A','log_PE-A')
#' @param sd_threshold Vector of standard deviation thresholds corresponding to the specified channels. Defaults to c(3, 3, 3, 3)
#' @param timebins Vector of the corresponding timebins for each fluorophore
#' @return A data frame of cell metadata with assignment columns (1 or -1 for each fluorophore/timebin) and a final timebin assignment column named 'group' 
#' @export
#' @examples
#' # Assuming 'my_data' is your flow cytometry data frame
#' result <- FACS_model(my_data)
#' result <- FACS_model(my_data, fluorophores = c('log_BV711-A','log_BUV737-A','log_BB515-A','log_PE-A'), sd_threshold =  c(3,3,3,3),
#'                      timepoints = c("12H","24H","36H","48H"))
#' @section Timebin Assignment
#' @import stats
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import MASS
#' @import dplyr

FACS_model <- function(well_fcs_mc, fluorophores = c('log_BV711.A','log_BUV737.A','log_BB515.A','log_PE.A'), sd_threshold = c(3,3,3,3),
                      timebins = c("12H","24H","36H","48H")) {
  
  channels = fluorophores
  # Initialize data.frame to store time_groups
  time_groups = data.frame(row.names = row.names(well_fcs_mc))

  # Filter for unstained samples
  unstained = well_fcs_mc$Stain == 'nonstained'
  
  # Loop over each channel to fit the model and classify cells
  for (i in seq_along(channels)){
    #frm <- as.formula(paste0("`",channels[i],"` ~ (`FSC.A` + `FSC.W` + `FSC.H` + `SSC.A` + `SSC.W` + `SSC.H` + `log_APC.A`)^2"))
    frm <- as.formula(paste0("`",channels[i],"` ~ (`FSC.A` + `FSC.W` + `FSC.H` + `SSC.A` + `SSC.W` + `SSC.H` + `log_APC.A`)^2"))
    g_model <- glm(frm,
                   family = gaussian(link = "identity"),
                   data = well_fcs_mc[unstained,])

    fit_model <- MASS::fitdistr((predict(g_model, type = "response") - well_fcs_mc[unstained, channels[i]]), "normal")
    bound = fit_model$estimate[1] - sd_threshold[i] * fit_model$estimate[2]

    # Calculate residuals
    x = predict(g_model, well_fcs_mc, type = "response")
    y = well_fcs_mc[, channels[i]]
    res = x - y
    res[is.na(res)] = 0

    # Classify cells based on residual values
    time_groups[, channels[i]] = 0
    time_groups[res < bound, channels[i]] = 1
    time_groups[res > bound, channels[i]] = -1
  }
  
  # Create groups based on channel conditions
  time_groups$group = 'Nan'
  time_groups$group[Reduce(`&`, lapply(channels, function(ch) time_groups[, ch] == -1))] = 'Negative'
  time_groups$group[time_groups[, channels[4]] == 1] = timebins[4]#'48H'
  time_groups$group[time_groups[, channels[3]] == 1] = timebins[3]#'36H'
  time_groups$group[time_groups[, channels[2]] == 1] = timebins[2]#'24H'
  time_groups$group[time_groups[, channels[1]] == 1] = timebins[1]#'12H'

  colnames(time_groups) = c(timebins, "group")
  well_fcs_mc = cbind(well_fcs_mc, time_groups)#$group
  return(well_fcs_mc)
}
                                       
#' plot_FACS_eval: Plot FACS data after assignment and check if assigned variable are good
#'
#' @param well_fcs_mc Data frame containing FACS data
#' @param fluorophores Vector of channels to be plotted
#' @param timebins Vector of the corresponding timebins for each fluorophore
#' @param y_axis String of fluorophore to be plotted against in the y-axis. Defaults to 'log_APC-A'
#' @return A ggplot object displaying the data
#' @export
#' @examples
#' # Assuming 'data' is your FACS data and 'assignment_data' contains the assignments
#' plot_FACS_eval_3(data, assignment_data)
#' plot_FACS_eval_3(data, assignment_data, channels = c('custom_channel1', 'custom_channel2'))
plot_FACS_eval <- function(well_fcs_mc, fluorophores = c('log_BV711.A','log_BUV737.A','log_BB515.A','log_PE.A'),
                             timebins =  c("12H","24H","36H","48H"), y_axis = 'log_APC.A') {
  df = well_fcs_mc
  channels = fluorophores
  plots = list()
  df$Stain = rep("stained", nrow(df))
  df$Stain[which(df$group=="Negative")] = "unstained"
  for (i in seq_along(channels)) {
    
    channel = channels[i]
    print(channel)
    timebin = timebins[i]
    g = ggplot(df, aes_string(x = channel, y= y_axis, color = "Stain")) + #x = paste0('`', channel, '`'), y = paste0('`', y_axis, '`')
        geom_point(size = .5, alpha = .1) +
        ggtitle(channel) +
        zmanr_theme
    
    g = g + geom_density_2d(data = df, aes_string(x = channel, y= y_axis, color = factor(df[, timebin])), alpha = .25) +
        scale_color_manual(values = c("black", "darkblue", "#008080", "#FF6F61"))
    
    plots[[i]] = g
  }
  
  options(repr.plot.width = 15, repr.plot.height = 10)
  
  ggarrange(plotlist = plots, ncol = 2, nrow = 2)
  
}