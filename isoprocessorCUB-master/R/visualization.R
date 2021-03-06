
#' Visualize reference peaks
#'
#' @note: refine and write tests
#' @param group_id first column in group id is used for x axis
#' @param within_group whether to visualize the deviation within the specified group or across all groups
#' @export
iso_visualize_ref_peaks <- function(dt, is_ref_condition, ratio = default(ratio), group_id = default(group_id), within_group = FALSE) {

  if (missing(dt)) stop("no data table supplied", call. = FALSE)
  if (missing(is_ref_condition)) stop("no condition as to what constitutes a reference peak provided", call. = FALSE)
  dt_cols <- get_column_names(!!enquo(dt), ratio = enquo(ratio), group_id = enquo(group_id), n_reqs = list(group_id = "+"))
  ref_quo <- enquo(is_ref_condition)

  refs <- dt %>%
    filter(!!ref_quo) %>%
    rename(ratio = !!as.name(dt_cols$ratio)) %>%
    mutate(delta_deviation = (ratio/mean(ratio) - 1) * 1000) %>%
    # calculate deviation from interrun ratio
    group_by(!!!cols_to_quos(dt_cols$group_id)) %>%
    mutate(
      run_delta_deviation = (ratio/mean(ratio) - 1) * 1000,
      ref_peak_nr = 1:n()) %>% ungroup() %>%
    mutate(ref_peak_nr = factor(ref_peak_nr))

  # plot reference peaks relative to total average
  p <- refs %>%
    ggplot() +
    aes_string(dt_cols$group_id[1], "delta_deviation", fill = "ref_peak_nr") +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "Analysis", fill = "Reference peak", y = "Deviation from total average [permil]")

  if (!within_group) {
    # across all samples
    return(p)
  } else {
    # inter-group average
    p %+%
      aes(y = run_delta_deviation) +
      labs(y = "Deviation within group [permil]")
  }
}

#' Visualize results for delta calibration
#' @inheritParams iso_print_data_table
#' @param x the x axis for the calibration parameters, by default a datetime but can also work as text
#' @param include_from_summary which parameters from the fit summary to include
#' @param date_breaks what breaks to use for the x axis if it is a datetime
#' @param date_labels datetime label pattern for x axis if it is a datetime
#' @export
iso_visualize_delta_calib_fits <- function(dt, x = default(file_datetime), color = NULL, shape = NULL, size = NULL,
                                           date_breaks = "2 hours", date_labels = "%d %b %H:%M",
                                           include_from_summary = c(adj.r.squared, deviance)) {

  # safety checks
  if (missing(dt)) stop("no data table supplied", call. = FALSE)
  if (length(missing <- setdiff(c("calib_delta_coefs", "calib_delta_summary"), names(dt))) > 0) {
    glue("missing columns in data table: '{collapse(missing, \"', '\")}'. Make sure to run iso_calibrate_delta() first and keep these columns.") %>%
      stop(call. = FALSE)
  }

  # pull out all coefficients (all for now, should always show all)
  calib_coefs <- dt %>% select(-calib_delta_summary) %>%
    iso_unnest_delta_calib_coefs(select = c(term, estimate, std.error))

  # pull out summary
  select_quo <- enquo(include_from_summary)
  calib_summary <- dt %>% select(-calib_delta_fit, -calib_delta_coefs) %>%
    iso_unnest_delta_calib_summary(select = !!select_quo)
  cs_cols <- get_column_names(calib_summary, select = select_quo, n_reqs = list(select = "*"))

  # check if any summary should be included
  if (length(cs_cols$select) > 0) {
    calib_summary <- gather(calib_summary, term, estimate, !!!cols_to_quos(cs_cols$select))
    visualization_data <- bind_rows(calib_coefs, calib_summary)
  } else {
    visualization_data <- calib_coefs
  }

  # check column availability
  vis_cols <- get_column_names(visualization_data, x = enquo(x), color = enquo(color), shape = enquo(shape), size = enquo(size),
                              n_reqs = list(color = "?", shape = "?", size = "?"))

  # make sure we only plot unique data sets
  visualization_data <- unique(visualization_data)

  # generate plot
  p <- visualization_data %>%
    # make sure term factor is in order (i.e. paneling)
    mutate(term = as_factor(term)) %>%
    filter(!is.na(estimate)) %>%
    ggplot() +
    aes_q(x = as.name(vis_cols$x), y = as.name("estimate"),
          color = if(!is_empty(vis_cols$color)) as.name(vis_cols$color) else NULL,
          shape = if(!is_empty(vis_cols$shape)) as.name(vis_cols$shape) else NULL,
          size = if(!is_empty(vis_cols$size)) as.name(vis_cols$size) else NULL
        ) +
    geom_errorbar(data = function(df) filter(df,!is.na(std.error)), map = aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0, size = 1) +
    { if (is_empty(vis_cols$size)) geom_point(size = 4) else geom_point() } +
    facet_wrap(~term, ncol = 1, scales = "free_y") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(plot.margin = margin(10, 5, 15, 20)) +
    labs(y="")

  # type of x axis
  if (is(visualization_data[[vis_cols$x]], "POSIXct")) {
    p <- p +
      scale_x_datetime(date_breaks = date_breaks, date_labels = date_labels) +
      labs(x = "")
  }

  return(p)
}

#' Visualize the data
#'
#' General purpose visualization function.
#'
#' @inheritParams iso_map_peaks
#' @inheritParams iso_visualize_delta_calib_fits
#' @param y which columns to visualize
#' @param group what to group by, multiple columns allowed (combine with c(...)), usually not necessary if color or other grouping is defined
#' @param lines whether to plot lines
#' @param points whether to plot points
#' @export
#' @note should probably make sure that the default columns for gather 'panel' and 'value' do not exist...
#' @note it would be great to allow renaming of the columns via this (especially the y column)
iso_visualize_data <- function(dt, x, y, group = NULL, color = NULL, shape = NULL, size = NULL, linetype = NULL, lines = TRUE, points = FALSE) {
  # safety checks
  if (missing(dt)) stop("no data table supplied", call. = FALSE)
  if (missing(x)) stop("have to provide an x to plot", call. = FALSE)
  if (missing(y)) stop("have to provide at least one column to plot", call. = FALSE)
  vis_cols <- get_column_names(!!enquo(dt), x = enquo(x), y = enquo(y), group = enquo(group),
                              color = enquo(color), shape = enquo(shape), linetype = enquo(linetype), size = enquo(size),
                              n_reqs = list(y = "+", group = "*", color = "?", shape = "?", linetype = "?", size = "?"))

  # gather data
  visualization_data <- gather(dt, panel, value, !!!cols_to_quos(vis_cols$y))

  # group quos (cols_to_quos does not work in this instance since it )
  group_quos <- quo(str_c(!!!cols_to_symbol_list(vis_cols$group)))

  # generate plot
  p <- visualization_data %>%
    # make sure panel factor is in order
    mutate(panel = as_factor(panel)) %>%
    filter(!is.na(value)) %>%
    ggplot() +
    aes_q(x = as.name(vis_cols$x), y = as.name("value"),
          group = if(!is_empty(vis_cols$group)) group_quos else NULL,
          color = if(!is_empty(vis_cols$color)) as.name(vis_cols$color) else NULL,
          shape = if(!is_empty(vis_cols$shape)) as.name(vis_cols$shape) else NULL,
          linetype = if(!is_empty(vis_cols$linetype)) as.name(vis_cols$linetype) else NULL,
          size = if(!is_empty(vis_cols$size)) as.name(vis_cols$size) else NULL
    ) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(plot.margin = margin(10, 5, 15, 20)) +
    labs(y="")

  # lines and points
  if (lines) {
    p <- p + geom_line()
  }
  if (points) {
    p <- p + { if (is_empty(vis_cols$size)) geom_point(size = 4) else geom_point() }
  }

  return(p)
}
