#' Conduct HYSPLIT trajectory runs
#'
#' The function executes single/multiple forward or backward HYSPLIT trajectory
#' runs using specified meteorological datasets.
#'
#' @param lat,lon,height The receptor position in terms of latitude and
#'   longitude (both in decimal degrees), and height in meters above ground
#'   level.
#' @param duration The duration of each model run (whether it is in the forward
#'   direction or running backwards) in hours.
#' @param days A vector of days that the model will run. This is combined with
#'   the `daily_hours` to produce a series of date-times.
#' @param daily_hours A vector of daily hours for initiations of runs across the
#'   given `days`. Use values from from `0` to `23`.
#' @param direction An option to select whether to conduct the model in the
#'   `"forward"` (default) or `"backward"` directions.
#' @param met_type The type of meteorological data files to use. The options
#'   are: `"reanalysis"` (NCAR/NCEP global reanalysis data, the default),
#'   `"gdas1"` and `"gdas0.5"` (Global Data Assimilation System 1-degree and
#'   0.5-degree resolution data), `"narr"` (North American Regional Reanalysis),
#'   `"gfs0.25"` (Global Forecast System 0.25 degree data), and `"nam12"` (North
#'   American Mesoscale Forecast System, 12-km/6-hour resolution data).
#' @param vert_motion A numbered option to select the method used to simulation
#'   vertical motion. The methods are: `0` (input model data), `1` (isobaric),
#'   `2` (isentropic), `3` (constant density), `4` (isosigma), `5` (from
#'   divergence), `6` (remap MSL to AGL), `7` (average data), and `8` (damped
#'   magnitude).
#' @param model_height The upper limit of the model domain in meters.
#' @param extended_met An option to report additional meteorological data along
#'   each output trajectory.
#' @param config A configuration list serves to internally generate the
#'   `SETUP.CFG` file. This list can be easily created by using the
#'   [set_config()] function. If `NULL`, then the default configuration list
#'   will be generated.
#' @param ascdata An ascdata list that will be used to create the `ASCDATA.CFG`
#'   file. This list can be provided through use of the [set_ascdata()]
#'   function. If `NULL`, then the default ascdata list will be generated.
#' @param traj_name An optional, descriptive name for the output file
#'   collection.
#' @param binary_path An optional path to a HYSPLIT trajectory model binary.
#'   When not specified, the model binary will be chosen from several available
#'   in the package (based on the user's platform).
#' @param met_dir An optional file path for storage and access of meteorological
#'   data files.
#' @param exec_dir An optional file path for the working directory of the model
#'   input and output files.
#' @param clean_up An option to make the `exec_dir` directory clean after
#'   completion of all trajectory runs. By default, this is set to `TRUE`.
#'   
#' @examples
#' \dontrun{
#' library(lubridate)
#' 
#' # Run a trajectory model 4 times a day
#' # for 6 days in 2012 using NCEP/NCAR
#' # reanalysis data
#' trajectory <-
#'   hysplit_trajectory(
#'     lat = 50.108,
#'     lon = -122.942,
#'     height = 100,
#'     duration = 48,
#'     days = seq(
#'       lubridate::ymd("2012-02-22"),
#'       lubridate::ymd("2012-02-27"),
#'       by = "1 day"
#'     ),
#'     daily_hours = c(0, 6, 12, 18)
#'   )
#' }
#' 
#' @export
hysplit_trajectory <- function(geo_df,
                               duration = 0,
                               date = NULL,
                               hour = 0,
                               direction = "forward",
                               met_type = "reanalysis",
                               vert_motion = 0,
                               model_height = 20000,
                               extended_met = FALSE,
                               temp_only = FALSE,
                               config = NULL,
                               ascdata = NULL,
                               traj_name = NULL,
                               binary_path = NULL,
                               met_dir = NULL,
                               exec_dir = NULL,
                               softrun = NULL, # This is not evaluated, yet
                               clean_up = TRUE) {
  daily_hours <- hour
  # geo_df <- geo_df %>%
  #   select(hour, date) %>%
  #   unique() %>%
  #   mutate(receptor = row_number()) %>%
  #   merge(run_df)
  
  ## TEMPORARY
  met_suff <- ifelse(met_type == 'nams', '_hysplit.t00z.namsa', '_nam12')
  
  # If the execution dir isn't specified, use the working directory
  if (is.null(exec_dir)) exec_dir <- getwd()
  else {
    if (!dir.exists(exec_dir)) {dir.create(exec_dir)}
  }
  
  # If the meteorology dir isn't specified, use the working directory
  if (is.null(met_dir)) met_dir <- paste0(getwd(), '/meteorology')
  else {
    if (!dir.exists(met_dir)) {dir.create(met_dir)}
  }
  
  #if (!is.null(geo_df)) {days <- unique(geo_df$date)}
  
  days <- as.Date(date)
  
  all_met_files <- list.files(met_dir)
  met_file_check <- get_daily_filenames(
    days, duration, direction, suffix = met_suff
  )
  
  infolder <- met_file_check %in% all_met_files
  if (!all(infolder)) {
    w <- which(!(infolder))
    stop('Missing the following met files: \n\n',
         paste(met_file_check[w], collapse = '\n'),
         '\n\nUse download_met_files()')
  }
  
  # Set the path for the `hyts_std` binary file
  binary_path <- 
    set_binary_path(
      binary_path = binary_path,
      binary_name = "hyts_std"
    )
  
  # Get the system type
  system_type <- get_os()
  
  # Generate name of output folder
  if (is.null(traj_name)) {
    folder_name <- paste0("traj-", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
  } else if (!is.null(traj_name)) {
    folder_name <- traj_name
  }
  
  if (is.null(config)) {
    config_list <- set_config()
  } else {
    config_list <- config
  }
  
  if (is.null(ascdata)) {
    ascdata_list <- set_ascdata()
  } else {
    ascdata_list <- ascdata
  }
  
  # Modify the default `SETUP.CFG` file when the option for extended
  # meteorology is `TRUE`
  if (isTRUE(extended_met)) {
    
    tm_names <-
      config_list %>%
      names() %>%
      vapply(
        FUN.VALUE = logical(1),
        USE.NAMES = FALSE,
        FUN = function(y) y %>% tidy_grepl(ifelse(temp_only, 
                                                  "^tm_tamb",
                                                  "^tm_"))
      ) %>%
      which()
    
    config_list[tm_names] <- 1
  }
  
  # Write the config and ascdata lists to files in
  # the `exec` directory
  config_list %>% write_config_list(dir = exec_dir)
  ascdata_list %>% write_ascdata_list(dir = exec_dir)
  
  # Stop function if there are vectors of different
  # length for `lat` and `lon`
  # if (length(lat) != length(lon)) {
  #   stop("The coordinate vectors are not the same length.", call. = FALSE)
  # }
  
  # Download any necessary meteorological data files
  # and return a vector of the all files required
  met_files <- 
    download_met_files(
      met_type = met_type,
      days = days,
      duration = duration,
      direction = direction,
      met_dir = met_dir
    )
  
  # Get vector of receptor indices
  
  start_year_GMT <- to_short_year(date)
  start_month_GMT <- to_short_month(date)
  start_day_GMT <- to_short_day(date)
  
  if (inherits(daily_hours, "numeric")) {
    daily_hours <- formatC(sort(daily_hours), width = 2, flag = 0)
  }
  
  output_filename <-
    get_traj_output_filename(
      traj_name = traj_name,
      direction = direction,
      year = start_year_GMT,
      month = start_month_GMT,
      day = start_day_GMT,
      hour = daily_hours,
      duration = duration
    )
  
  write_traj_control_file(
    start_year_GMT = start_year_GMT,
    start_month_GMT = start_month_GMT,
    start_day_GMT = start_day_GMT,
    start_hour_GMT = daily_hours,
    geo_df = geo_df,
    direction = direction,
    duration = duration,
    vert_motion = vert_motion,
    model_height = model_height,
    met_files = met_files,
    output_filename = output_filename,
    system_type = system_type,
    met_dir = met_dir,
    exec_dir = exec_dir
  )
}
