
## Make sure the "days" variable is in Date format
#' Download Meteorology Files
#'
#' @param days Dates for which to pull meteorology data
#' @param duration Duration of desired HYSPLIT run; determines how many dates
#' before/after "days" are pulled
#' @param direction Direction of desired HYSPLIT run; determines whether to pull
#' additional days before or after "days"
#' @param met_dir Path to directory where meteorology files are to be written
#'
#' @return Returns names of met files; downloads necessary meteorology files.
#' @export
#'
#' @examples
download_met_files <- function(days,
                               duration,
                               direction,
                               met_dir = paste0(getwd(), '/meteorology')) {
  options(timeout = 800)
  
  if (!dir.exists(met_dir)) {
    dir.create(met_dir)
  }
  
  metfiles <- paste0(met_dir, '/', list.files(met_dir, pattern = 't00z.namsa'))
  metfile_sizes <- file.info(metfiles)$size
  true_sizes <- sort(unique(metfile_sizes), decreasing = TRUE)[1:2]
  min_size <- min(true_sizes)
  small_file_inds <- which(metfile_sizes < min_size)
  #file.remove(metfiles[small_file_inds])
  if (length(small_file_inds) > 0) {
    warning(paste0('Meteorology directory contains ',
                   length(small_file_inds),
                   ' met files which are smaller than ',
                   min_size/1000,
                   'KB. PLEASE check that all met files are complete before proceeding!!'))
  }
  
  get_daily_filenames(
    days = days,
    duration = duration,
    direction = direction,
    suffix = "_hysplit.t00z.namsa"
  ) %>%
    get_met_files(
      path_met_files = met_dir,
      ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/nams"
    )
  
}