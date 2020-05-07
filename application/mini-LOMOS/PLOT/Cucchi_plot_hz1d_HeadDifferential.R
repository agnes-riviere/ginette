

#' Plots temperature measurement timeseries
#'
#' @param param the vector of parameters
#' @param df_hz1d_t the dataframe of measurements
#' @return a ggplot object
#' @export
plot_hz1d_t <- function(param,df_hz1d_t){

  df_z_idx <- HZinv::get_level_z_fromParam(param = param)

  # create named vector for legend labels
  vect_labels <- as.character(df_z_idx$z)
  names(vect_labels) <- as.character.factor(df_z_idx$z_idx)

  # final plot
  g_res <-
    ggplot() +
    geom_line(data = df_hz1d_t,
              mapping = aes(x=t_time,y=temperature,color=z_idx)) +
    scale_color_manual(
      values=col_z, # internal value for col_z defined in data-raw/internal_colors.R
      labels = vect_labels) +
    labs(x="",y="T (in C)", color="depth\n(in m)") +
    theme_bw()

  return(g_res)

}

#' Plots head differential and temperature measurement timeseries from dataframe
#' from csv file
#'
#' @param df_measGinette the dataframe of measurements (as outputted by
#'   read_measGinette)
#' @param dates_range range of dates for x-axis
#' @return a ggplot object
#' @export
plot_fromMeas <- function(df_measGinette,dates_range){

  g_meas <-
    cowplot::plot_grid(
      HZinv::plot_head_differential_fromMeas(df_measGinette = df_measGinette,
                                             dates_range = dates_range),
      HZinv::plot_temperature_fromMeas(df_measGinette = df_measGinette,
                                       dates_range = dates_range),
      nrow = 2,rel_heights = c(4,5))

  return(g_meas)

}

#' Plots temperature measurement timeseries from dataframe from csv file
#'
#' @param df_measGinette the dataframe of measurements
#'  (as outputted by read_measGinette)
#' @param dates_range range of dates for x-axis
#' @return a ggplot object
#' @export
plot_temperature_fromMeas <- function(df_measGinette,dates_range=NULL){

  # define vector of labels
  df_labels <-
    unique(dplyr::left_join(x = data.frame(variable=names(HZinv:::col_z)),
                            y = df_measGinette[,c('variable','depth_m')]))
  vect_labels <- df_labels$depth_m
  names(vect_labels) <- df_labels$variable

  if(is.null(dates_range)){
    dates_range = range(df_measGinette$t_time,na.rm = T)
  }

  g_meas_temperature <-
    ggplot(data = df_measGinette) +
    geom_line(data = subset(x = df_measGinette,
                            subset = type == "temperature_C"),
              mapping = aes(x=t_time,y=value,
                            col=factor(variable,levels = names(HZinv:::col_z))))  +
    labs(x="",y="T (in C)", color="depth below streambed (in m)") +
    scale_x_datetime(labels=date_format("%d %b %Y"),
                     limits=dates_range) +
    scale_color_manual(
      values = HZinv:::col_z, # internal value for col_z defined in data-raw/internal_colors.R
      labels = vect_labels) +
    theme_bw() +
    theme(legend.position = "top")

  return(g_meas_temperature)

}

#' Plots head differential measurement timeseries from dataframe from csv file
#'
#' @param df_measGinette the dataframe of measurements
#'  (as outputted by read_measGinette)
#' @param dates_range range of dates for x-axis
#' @return a ggplot object
#' @export
plot_head_differential_fromMeas <- function(df_measGinette,dates_range = NULL){

  if(is.null(dates_range)){
    dates_range = range(df_measGinette$t_time,na.rm = T)
  }

  g_meas_head_differential <-
    ggplot(data = subset(x = df_measGinette,
                         subset = variable == "head_differential_m")) +
    geom_line(mapping = aes(x = t_time,y = value)) +
    expand_limits(y = 0) +
    geom_hline(mapping = aes(yintercept = 0),linetype = "dashed") +
    labs(x="",y = expression(Delta*'H = H'['HZ'] *'- H'['riv'] * ' (in m)')) +
    scale_x_datetime(labels=date_format("%d %b %Y"),
                     limits=dates_range) +
    theme_bw()

  return(g_meas_head_differential)

}




