met_file_check <- function(met_dir, dates, duration, direction) {
  all_met_files <- list.files(met_dir, pattern = "_hysplit.t00z.namsa")
  met_file <- get_daily_filenames(as.Date(dates), duration, direction, 
                                  suffix = "_hysplit.t00z.namsa")
  
  infolder <- met_file %in% all_met_files
  if (!all(infolder)) {
    w <- which(!(infolder))
    stop('Missing the following met files: \n\n',
         paste(met_file[w], collapse = '\n'),
         '\n\nUse sbwsplitr::download_met_files()')
  }
  return(met_file)
}