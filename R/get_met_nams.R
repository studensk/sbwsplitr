#' Get NAMS meteorology data files
#'
#' Downloads NAMS meteorology data files from the NOAA FTP server and saves
#' them to a specified folder. Files can be downloaded by specifying a list of
#' filenames (in the form of `"{YYYY}{MM}{DD}_nam12"`).
#' 
#' @inheritParams get_met_gdas1
#' 
#' @export
get_met_nams <- function(days,
                         duration,
                         direction,
                         path_met_files) {
  
  get_daily_filenames(
    days = days,
    duration = duration,
    direction = direction,
    suffix = "_hysplit.t00z.namsa"
  ) %>%
    get_met_files(
      path_met_files = path_met_files,
      ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/nams"
    )
}
