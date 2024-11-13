#' Conduct multiple HYSPLIT trajectory runs
#'
#' The function executes multiple forward or backward HYSPLIT trajectory
#' runs using specified meteorological datasets.
#'
#' @param run_df A data frame whose rows contain the information for each 
#'   trajectory run. This data frame should contain the following columns:
#'   "lat", "lon", "height", "date" and "hour"
#' @param days A vector of days that the model will run. This is combined with
#'   the `hours` to produce a series of date-times.
#' @param hours A vector of daily hours for initiations of runs across the
#'   given `days`. Use values from from `0` to `23`.
#' @param lat,lon,height The receptor position in terms of latitude and
#'   longitude (both in decimal degrees), and height in meters above ground
#'   level.
#' @param duration The duration of each model run (whether it is in the forward
#'   direction or running backwards) in hours.
#' @param direction An option to select whether to conduct the model in the
#'   `"forward"` (default) or `"backward"` directions.
#' @param model_height The upper limit of the model domain in meters.
#' @param vert_motion A numbered option to select the method used to simulation
#'   vertical motion. The methods are: `0` (input model data), `1` (isobaric),
#'   `2` (isentropic), `3` (constant density), `4` (isosigma), `5` (from
#'   divergence), `6` (remap MSL to AGL), `7` (average data), and `8` (damped
#'   magnitude).
#' @param met_dir An optional file path for storage and access of meteorological
#'   data files. 
#' @param exec_dir An optional file path for the working directory of the model
#'   input and output files.
#' @param binary_path An optional path to a HYSPLIT trajectory model binary.
#'   When not specified, the model binary will be chosen from several available
#'   in the package (based on the user's platform).
#' @param csv_folder An optional, descriptive name for the output file
#'   collection.
#' @param clean_up An option to make the `exec_dir` directory clean after
#'   completion of all trajectory runs. By default, this is set to `TRUE`.
#' @param rdf_write A character string indicating the name of the output file 
#'   for the input data. If NULL, no file is written.
#' @param traj_write A character string indicating the name of the file 
#'   for the trajectory output. If NULL, no file is written.
#' @param plot TRUE or FALSE; produce trajectory plots?
#' @param vbug Moth velocity in \eqn{m/s^2}
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
multiple_trajectories <- function(run_df = NULL,
                                  days = NULL,
                                  hours = NULL,
                                  lat = NULL,
                                  lon = NULL,
                                  height = NULL,
                                  duration = 9,
                                  direction = 'forward',
                                  model_height = 20000,
                                  vert_motion = 0,
                                  met_dir = getwd(),
                                  exec_dir = getwd(), 
                                  binary_path = NULL,
                                  csv_folder = 'traj_output',
                                  clean_up = TRUE,
                                  rdf_write = NULL, 
                                  traj_write = NULL,
                                  plot = FALSE, 
                                  vbug = 2.5, ...) {
  
  if (is.null(run_df)) {
    run_df <- make_run_df(lat, lon, height, days, hours)
  }
  
  # If the execution dir isn't specified, use the working directory
  if (!dir.exists(exec_dir)) {dir.create(exec_dir)}
  
  output_path <- file.path(exec_dir, csv_folder)
  if (!dir.exists(output_path)) {dir.create(output_path)}
  
  # If the meteorology dir isn't specified, use the working directory
  #if (!dir.exists(met_dir)) {dir.create(met_dir)}
  
  if (!is.null(run_df)) {days <- unique(run_df$date)}
  
  # days <- sort(as.Date(days))
  days <- as.character(sort(days))
  met_file <- met_file_check(met_dir, days, duration, direction)
  
  if (is.null(binary_path)) {
    binary_path_set <-
      set_binary_path(
        binary_path = binary_path,
        binary_name = "hyts_std"
      )
  }
  
  run_tbl <- run_df %>%
    dplyr::select(lat, lon) %>%
    dplyr::distinct() %>%
    dplyr::mutate(site = dplyr::row_number()) %>%
    merge(run_df) %>%
    dplyr::arrange(date, site) %>%
    dplyr::mutate(traj_id_full = paste('traj', date, height,
                                       site, hour, sep = '_'),
                  traj_id = dplyr::row_number())
  
  if (!is.null(rdf_write)) {
    write.csv(run_tbl, 
              paste0(file.path(exec_dir, rdf_write), '.csv'),
              row.names = FALSE)
  }
  for (d in days) {
    date_df <- subset(run_tbl, date == d)
    sites <- unique(date_df$site)
    traj.lst <- lapply(sites, function(s) {
      ds_df <- subset(date_df, site == s)
      run_dir <- file.path(exec_dir, paste(d, s, sep = '_'))
      dir.create(run_dir)
      for (r in 1:nrow(ds_df)) {
        run_vals <- as.list(ds_df[r,])
        traj <- hysplit_trajectory(run_vals = run_vals,
                                   duration = duration,
                                   direction = direction,
                                   model_height = model_height,
                                   vert_motion = vert_motion,
                                   vbug = vbug,
                                   met_dir = met_dir,
                                   directory = run_dir, 
                                   bin_path = binary_path_set,
                                   return_traj = FALSE)
      }
      
      # all.files <- list.files(receptor_dir)
      # traj.filename <- all.files[grep('traj', all.files)]
      # traj.file.orig <- file.path(receptor_dir, traj.filename)
      # traj.file <- file.path(output_path, traj.filename)
      # file.copy(traj.file.orig, traj.file)
      
      
      traj_tbl <-
        trajectory_read(output_folder = run_dir) %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(
          lat_i = lat_i,
          lon_i = lon_i,
          height_i = height_i,
          date_i = as.Date(traj_dt_i),
          site_i = s,
          traj_id_full = paste('traj', date_i, height_i,
                               site_i, hour_i, sep = '_')
        )
      write.file <- paste0('traj_', as.character(d), '_', s, '.csv')
      write.path <- file.path(output_path, write.file)
      
      write.csv(traj_tbl, write.path, row.names = FALSE)
      unlink(run_dir, recursive = TRUE, force = TRUE)
      
    })
  }
  
  all.traj.files <- list.files(output_path)
  traj.lst <- lapply(all.traj.files, function(file) {
    read.csv(file.path(output_path, file))
  }) 
  traj.df <- bind_rows(traj.lst) 
  run.df.mg <- run_tbl %>%
    select(traj_id_full, traj_id) %>%
    merge(traj.df)
  
  if (!is.null(traj_write)) {
    write.csv(run.df.mg, file.path(exec_dir, paste0(traj_write_name, '.csv')), 
              row.names = FALSE)
  }
  
  if (clean_up) {unlink(output_path, recursive = TRUE, force = TRUE)}
  
  if (plot) {
    print(plot_trajectories(run.df.mg, ...))
  }
  return(run.df.mg)
}
