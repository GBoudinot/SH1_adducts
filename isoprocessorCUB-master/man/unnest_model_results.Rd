% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{unnest_model_results}
\alias{unnest_model_results}
\title{Unnest model results}
\usage{
unnest_model_results(dt, model_results, select = everything(),
  keep_remaining_nested_data = FALSE, keep_other_list_data = FALSE)
}
\arguments{
\item{model_results}{name of the model results column to unnest}

\item{select}{which fields from the selected model result column to unnest}

\item{keep_remaining_nested_data}{whether to keep any remaining parts of the partially unnested data}

\item{keep_other_list_data}{keep other list data columns (e.g. other data or model columns)}
}
\description{
Convenience functions for unnesting model results
}
