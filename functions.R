#' map peaks
#' @param raw data frame with raw data (data table export)
#' @param metadata_file file path to the metadata.xslx file
map_peaks <- function(raw, metadata_file, quiet = FALSE) {
  stopifnot(file.exists(metadata_file))
  stopifnot(all(c("Start", "Rt", "End") %in% names(raw)))
  
metadata_maps <- read_excel(metadata_file , sheet = "maps")
  metadata_maps[names(metadata_maps) %in% c("", "``")] <- NULL
  stopifnot(all(c("compound", "is_ref_peak") %in% names(metadata_maps)))
  
metadata_files <- read_excel(metadata_file, sheet = "files")
  metadata_files[names(metadata_files) %in% c("", "``")] <- NULL
  stopifnot(all(c("File", "map", "process") %in% names(metadata_files)))
  
  # read peak maps
  maps <<- metadata_maps %>% 
    filter(!is.na(compound)) %>% # ignore empty rows
    bind_rows(data_frame(compound = NA, is_ref_peak = "unknown")) %>% unique() %>% # add NA compound to map undefined peaks
    gather(map, Rt_target, starts_with("Rt")) %>% # gather maps
    mutate(map = sub("Rt\\.", "", map)) %>% # map name
    right_join(metadata_files %>% mutate(map = as.character(map)), by = "map") %>% # merge in file metadata
    filter(!(is.na(Rt_target) & !is.na(compound))) %>%  # remove if compound defined but not retention time
    filter(!is.na(map)) # ignore completely if map is not defined
  
  # combine peak map data and file data
  data_all <<- suppressMessages(
    raw %>% 
      # join in the specific peak maps
      left_join(maps, by = c("File")) %>% 
      # find the peak that the retention time window
      mutate(is_target = is.na(Rt_target) | (Start < Rt_target & End > Rt_target)) %>% 
      # figure out which peaks match multiple definitions and which definitions match multiple peaks
      group_by(File, Nr.) %>% mutate(n_matches = sum(is_target) - 1) %>% 
      group_by(File, Rt_target) %>% mutate(n_overlapping = ifelse(!is.na(compound), sum(is_target) - 1, 0)) %>% 
      ungroup() %>% 
      # filter out definitios that don't fit
      filter( (n_matches > 0 & is_target & !is.na(compound)) | (n_matches == 0 & is.na(compound))) %>% 
      # remove temp field is_target
      select(-is_target) %>% 
      # make sure listed peaks that have no matches are captured
      full_join(filter(maps, File %in% unique(raw$File))) %>% # dynamically join by all fields
      # arrange in useful order
      mutate(Rt_arrange = ifelse(!is.na(Rt_target), Rt_target, Rt)) %>% 
      arrange(File, Rt_arrange, Nr.) %>% select(-Rt_arrange) %>% 
      # remove redundant no peak placeholder
      filter(!(is.na(Nr.) & is.na(compound))) %>% 
      # remove all files that are not supposed to be processed
      filter(process == "yes") 
  )
  
  # only the assigned data
  data_matched <- data_all %>% filter(!is.na(Nr.), !is.na(compound))
  message(sprintf("INFO: %d peaks in %d files successfully assigned.", nrow(data_matched), length(unique(data_matched$File))))
  
  # mapping warnings
  if (!quiet) {
    if (nrow(map_missing <<- filter(data_all, is.na(map) | is.na(is_ref_peak)) %>% select(File) %>% distinct()) > 0) 
      message(sprintf("WARNING: %d file(s) have no (or undefined) peak map assignments\n%s", nrow(map_missing), paste(map_missing$File, collapse = "\n")))
    
    if (nrow(map_duplicates <<- filter(data_all, n_matches > 1 | n_overlapping > 0)) > 0) {
      message(sprintf("WARNING: %d peaks have duplicate assignments or are overlapping", nrow(map_duplicates)))
      print(map_duplicates)
    }
    
    if (nrow(map_unassigned <<- filter(data_all, !is.na(map), is.na(compound) | is.na(Nr.))) > 0) {
      message(sprintf("WARNING: %d peaks are unassigned or missing", nrow(map_unassigned)))
      print(map_unassigned)
    }
  }
  
  return(data_matched)
}