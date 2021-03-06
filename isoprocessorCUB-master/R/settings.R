#' @export
isoreader::iso_turn_info_messages_on

#' @export
isoreader::iso_turn_info_messages_off

# default parameters =====

# retrieve package settings, internal function, not exported
# first checks if the setting exists in the isoreader space, then in the isoprocessor
# @note: consider providing an option that indicates whether to return the quoted expression or a default value
# e.g. default(x) that returns quo(x) if x is not set vs. default(x, NULL) that returns quo(NULL) if x is not set
# @FIXME: alternatively set ALL available parameters in the initialize_options and force parameters to exist
default <- function(name) {
  if (missing(name)) return(quo())
  name_quo <- enquo(name)
  name <- name_quo %>% quos_to_text(variable = "setting")
  value <- isoreader:::default(!!as.name(name), allow_null = TRUE)
  if (is.null(value)) # not in isoreader settings
    value <- getOption(str_c("isoprocessor.", name))
  if (is.null(value)) { # not in normal isoprocessor settings
    func_params <- get_process_parameters()
    if (name %in% names(func_params))
      value <- func_params[[name]]
    else
      value <- name_quo
  }
  return(value)
}

# set package setting, internal function, not exported
# @note essetially same function as in isoreader except with the isoprocessor. prefix
set_default <- function(name, value, overwrite = TRUE) {
  if (overwrite || !str_c("isoprocessor.", name) %in% names(options()))
    options(list(value) %>% setNames(str_c("isoprocessor.", name)))
  return(invisible(value))
}

# processing default parameters =====

# get the default parameters
get_process_parameters <- function() {
  default("default_parameters")
}

#' Set default parameters
#' Define the default parameters for calculations used in this package.
#' @param data can be used to include this function call within a pipeline
#' @param ... named entries for the default parameters (standard or non-standard evaluation supported).
#' @return invisibly returns the passed in \code{data} object
#' Names are equivalent to function parameter names and
#' @family settings functions
#' @export
iso_set_default_process_parameters <- function(data = NULL, ...) {
  def_cols <- quos(...)
  existing_def_cols <- get_process_parameters()
  set_default("default_parameters", modifyList(existing_def_cols, def_cols))
  return(invisible(data))
}

#' @details \code{iso_reset_default_process_parameters} resets all default function parameters for this package.
#' @rdname iso_set_default_process_parameters
#' @export
iso_reset_default_process_parameters <- function(data = NULL) {
  initialize_options()
  return(invisible(data))
}

#' Get the current default processor parameters
#'
#' Retrieve a table with all default parameters for the isoprocessor package
#' (see \link[isoreader]{iso_get_default_reader_parameters} for the equivalent isoreader function).
#' To set processor parameters, see \code{\link{iso_set_default_processor_parameters}}.
#' To set messaging and caching parameters see \code{\link[isoreader]{iso_info_messages}}.
#' For a piping compatible version of this function, see \link{iso_show_default_processor_parameters}.
#' @inheritParams iso_set_default_process_parameters
#' @family settings functions
#' @export
iso_get_default_processor_parameters <- function(data = NULL) {
  regular_params <- data_frame(parameter = "quiet", value = as.character(default("quiet")))
  func_params <- get_process_parameters()
  if (length(func_params) > 0) {
    bind_rows(
      regular_params,
      data_frame(parameter = names(func_params), value = map_chr(func_params, quo_text))
    ) %>% return()
  } else
    return(regular_params)
}

#' Show the current default processor parameters
#'
#' Shows a table with the default function parameters for this package.
#' @inheritParams iso_set_default_process_parameters
#' @param func function to use for formatting the reader parameters table, e.g. \code{\link[knitr]{kable}}.
#' Note that if the output is in RMarkdown chunks, the chunk option must have \code{results="asis"} for the table to be correctly formatted.
#' @param ... additional parameters to forward to the \code{func} function
#' @family settings functions
#' @export
iso_show_default_processor_parameters <- function(data = NULL, print_func = default(print_func), ..., quiet = default(quiet)) {
  if (!quiet) message("Info: isoprocessor package current default parameters")

  print_func <- enquo(print_func) %>% resolve_defaults()
  if (!quo_is_null(print_func))
    print(do.call(eval_tidy(UQE(print_func)), args = c(list(x = iso_get_default_processor_parameters()), list(...))))
  else
    print(iso_get_default_processor_parameters())

  # for pipeline
  return(invisible(data))
}
