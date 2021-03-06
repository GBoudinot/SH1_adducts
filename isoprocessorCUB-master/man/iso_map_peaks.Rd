% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peak_mapping.R
\name{iso_map_peaks}
\alias{iso_map_peaks}
\title{Map peaks based on retention time}
\usage{
iso_map_peaks(dt, peak_maps, file_id = default(file_id),
  map_id = default(map_id), compound = default(compound),
  rt = default(rt), rt_start = default(rt_start),
  rt_end = default(rt_end), rt_prefix_divider = ":",
  is_identified = default(is_identified), is_missing = default(is_missing),
  is_ambiguous = default(is_ambiguous), n_matches = default(n_matches),
  n_overlapping = default(n_overlapping), quiet = default(quiet))
}
\arguments{
\item{dt}{data frame with peak data}

\item{peak_maps}{data frame with the peak map}

\item{file_id}{the column in dt that holds file id information}

\item{map_id}{the column in dt that indicates which map to use for which file}

\item{compound}{the column in peak_maps that holds compound information}

\item{rt}{the column in dt and colum prefix in peak_maps ("rt:...") that holds retention time information}

\item{rt_start}{the column in dt that holds start of peak retention times}

\item{rt_end}{the column in dt that holds end of peak retention times}

\item{rt_prefix_divider}{the divider after the retention time column prefix in peak_maps to identify the map id values (e.g. "rt:map_id_value")}

\item{is_identified}{the name for the new column which stores a logical TRUE/FALSE indicating peaks that have been successfully identified (includes missing peaks from the peak map!)
(note that this information could also be derived from !is.na(compound) but is provided for convenience)}

\item{is_missing}{the name for the new column which stores a logical TRUE/FALSE indicating peaks that are in the peak map definition but have no matching peak}

\item{is_ambiguous}{the name for the new column which stores a logical TRUE/FALSE indicating peaks that are ambiguous in their definition either because
they have multiple matches or because they overlap with other, overlapping peaks that were identified the same (note that this information could also be derived
from n_overlapping > 1 | n_matches > 1 but is provided for convenience)
(note that this information could also be derived from n_overlapping == 0 but is provided for convenience)}

\item{n_matches}{the name for the new column which stores the number of matches each peak has in the peak map}

\item{n_overlapping}{the name for the new column which stores the number of overlapping peaks that match the same peak definition}
}
\description{
This function makes it easy to map peaks based on peak maps. It reports all peaks including missing peaks and ambiguous peaks.
For routine downstream data processing, this function is usually followed by \code{filter(is_identified, !is_missing, !is_ambiguous)} to proceed only with
peaks that do not have any problems. Note that without this filter, one must proceed with great caution interpreting the ambiguous peaks.
}
\note{
: TODO - allow multiple file_id, not just one
}
