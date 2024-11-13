#' Conduct a single HYSPLIT trajectory run
#'
#' The function executes one forward or backward HYSPLIT trajectory
#' runs using specified meteorological datasets.
#'
#' @param run_vals List of starting values for the trajectory. Must contain lat,
#'   lon, height, date and hour.
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
hysplit_trajectory <- function(run_vals = NULL,
                               date = NULL,
                               hour = NULL,
                               lat = NULL,
                               lon = NULL,
                               height = NULL,
                               site = NULL,
                               duration = 9,
                               direction = 'forward',
                               model_height = 20000,
                               vert_motion = 0,
                               vbug = 2.5,
                               met_dir = file.path(getwd(), 'meteorology'),
                               directory = NULL, 
                               bin_path = NULL, 
                               return_traj = TRUE) {
  
  if (is.null(met_dir)) {stop('Please specify met_dir')}
  if (is.null(directory)) {stop('Please specify directory for HYSPLIT run')}
  if (is.null(bin_path)) {
    bin_path <-
      set_binary_path(
        binary_path = bin_path,
        binary_name = "hyts_std"
      )
  }
  
  if (!is.null(run_vals)) {
    lat <- run_vals$lat
    lon <- run_vals$lon
    height <- run_vals$height
    date <- run_vals$date
    hour <- run_vals$hour
    site <- run_vals$site
  }
  
  if (!dir.exists(directory)) {dir.create(directory)}
  
  else {
    null.val <- sapply(list(date, hour, lat, lon, height), is.null)
    if (any(null.val)) {
      w <- which(null.val)
      variables <- c('date', 'hour', 'lat', 'lon', 'height')
      stop(paste0('Value of ', variables[w[1]], ' is NULL'))
    }
  }
  
  met_file <- met_file_check(met_dir, date, duration, direction)
  
  system_type <- get_os()
  
  start_year_GMT <- to_short_year(date)
  start_month_GMT <- to_short_month(date)
  start_day_GMT <- to_short_day(date)
  
  # Sort daily starting hours if given as
  # numeric values
  if (inherits(hour, "numeric")) {
    hour <- formatC(sort(hour), width = 2, flag = 0)
  }
  
  start_hour_GMT <- hour
  full_year_GMT <- as.character(year(as.Date(date)))
  
  config_list <-  list(KMSL = 0,
                       tm_tpot = 1,
                       tm_tamb = 1,
                       tm_rain = 1,
                       tm_mixd = 1,
                       tm_relh = 1,
                       tm_terr = 1,
                       tm_dswf = 1,
                       vbug=vbug)
  
  ascdata_list <- set_ascdata()
  
  config_list %>% write_config_list(dir = file.path(directory))
  ascdata_list %>% write_ascdata_list(dir = file.path(directory))
  
  output_filename <-
    get_traj_output_filename(
      traj_name = NULL,
      site = site,
      direction = direction,
      year = start_year_GMT,
      month = start_month_GMT,
      day = start_day_GMT,
      hour = start_hour_GMT,
      lat = lat,
      lon = lon,
      height = height,
      duration = duration
    )
  
  
  # Write the CONTROL file
  write_traj_control_file(
    start_year_GMT = start_year_GMT,
    start_month_GMT = start_month_GMT,
    start_day_GMT = start_day_GMT,
    start_hour_GMT = start_hour_GMT,
    lat = lat,
    lon = lon,
    height = height,
    direction = direction,
    duration = duration,
    vert_motion = vert_motion,
    model_height = model_height,
    met_files = met_file,
    output_filename = output_filename,
    system_type = system_type,
    met_dir = met_dir,
    exec_dir = directory
  )
  
  # The CONTROL file is now complete and in the
  # working directory, so, execute the model run
  sys_cmd <-
    paste0(
      "(cd \"",
      directory,
      "\" && \"",
      bin_path,
      "\" ",
      to_null_dev(system_type = system_type),
      ")"
    )
  
  execute_on_system(sys_cmd, system_type = system_type)
  
  if (return_traj) {
    traj_tbl <-
      trajectory_read(output_folder = directory) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(
        lat_i = lat_i,
        lon_i = lon_i,
        height_i = height_i,
        date_i = as.Date(traj_dt_i),
        traj_id_full = paste('traj', date_i, height_i,
                             hour_i, sep = '_')
      )
    unlink(file.path(directory, output_filename), recursive = TRUE)
    return(traj_tbl)
  }
}