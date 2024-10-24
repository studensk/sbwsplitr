## Hysplit Runs ##
## run_df should have lat, lon, height, date and hour columns. If not in GMT,
##   there should also be a "timezone" column
##
## if not using run_df, it's assumed that all of the times are in GMT

#' Calculate HYSPLIT Trajectories
#'
#' @param run_df A data frame whose rows contain the information for each
#' trajectory run. This data frame should contain the following columns: "lat",
#' "lon", "height", "date" and "hour"
#' @param lat Latitude
#' @param lon Longitude
#' @param height Dispersal altitude
#' @param duration Number of hours for the trajectories to run
#' @param days Dates for the trajectory runs (as character or Date)
#' @param daily_hours Starting hours for the trajectory (numeric)
#' @param direction One of "forward" or "backward", denoting the direction of
#' the HYSPLIT run
#' @param vert_motion
#' @param model_height
#' @param extended_met TRUE or FALSE; should extra meteorological columns be
#' included in the output?
#' @param vbug Moth velocity
#' @param traj_name Name of output folder
#' @param binary_path1
#' @param met_dir Path to directory containing downloaded meteorological data
#' @param exec_dir Location of main directory
#' @param clean_up TRUE or FALSE; delete all run folders on completion?
#' @param local_time TRUE or FALSE; the input times are local relative to the
#' input latitude and longitude. If FALSE, it is assumed that times are in GMT.
#'
#' @return
#' @export
#'
#' @examples
hysplit_trajectory <- function(run_df = NULL,
                               lat = NULL,
                               lon = NULL,
                               height = 500,
                               duration = 9,
                               days = NULL,
                               daily_hours = 0,
                               direction = "forward",
                               vert_motion = 0,
                               model_height = 20000,
                               extended_met = TRUE,
                               vbug = 2.5,
                               traj_name = NULL,
                               binary_path1 = NULL,
                               met_dir = paste0(getwd(), '/meteorology'),
                               exec_dir = NULL,
                               clean_up = TRUE,
                               local_time = FALSE) {
  
  config <- list(KMSL = 0,
                 tm_tpot = 1,
                 tm_tamb = 1,
                 tm_rain = 1,
                 tm_mixd = 1,
                 tm_relh = 1,
                 tm_terr = 1,
                 tm_dswf = 1,
                 vbug=vbug)
  
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
  
  if (!is.null(run_df)) {days <- unique(run_df$date)}
  
  days <- as.Date(days)
  
  all_met_files <- list.files(met_dir)
  met_file_check <- get_daily_filenames(
    days, duration, direction, suffix = "_hysplit.t00z.namsa"
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
      binary_path = binary_path1,
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
  recep_file_path <- file.path(exec_dir, folder_name)
  
  if (is.null(config)) {
    config_list <- set_config()
  } else {
    config_list <- config
  }
  
  # if (is.null(ascdata)) {
  #   ascdata_list <- set_ascdata()
  # } else {
  #   ascdata_list <- ascdata
  # }
  
  ascdata_list <- set_ascdata()
  
  print('met cleared')
  
  ##  Generate a tibble of receptor sites
  # Stop function if there are vectors of different
  # length for `lat` and `lon`
  if (is.null(run_df)) {
    if (length(lat) != length(lon)) {
      stop("The coordinate vectors are not the same length.", call. = FALSE)
    }
    receptors_tbl <-
      dplyr::tibble(lat = lat, lon = lon) %>%
      dplyr::group_by(lat, lon) %>%
      tidyr::expand(height = height, hour = daily_hours, date = days) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(receptor = dplyr::row_number()) %>%
      dplyr::select(receptor, dplyr::everything())
  }
  
  else {
    receptors_tbl <- run_df %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(receptor = dplyr::row_number()) %>%
      dplyr::select(receptor, dplyr::everything())
    
  }
  
  if (local_time) {
    receptors_tbl$timezone <- tz_lookup_coords(lat = receptors_tbl$Latitude,
                                               lon = receptors_tbl$Longitude,
                                               method = 'accurate')
    dt <- as.POSIXct(paste0(receptors_tbl$date, '00:00'),
                     format = '%Y-%m-%d %H:%M',
                     tz = receptors_tbl$timezone)
    new_dt <- with(new, 'GMT')
  }
  
  # Get vector of receptor indices
  receptors <- seq(nrow(receptors_tbl))
  
  cores <- parallel::detectCores()
  max_clusters <- floor(cores*2/3)
  clusters <- pmin(max_clusters, max(receptors))
  print(paste0(clusters, ' clusters'))
  
  print('initialize clusters')
  cl <- makeCluster(clusters)
  clusterEvalQ(cl, {
    library(tidyverse)
  })
  clusterExport(cl, c("receptors_tbl", "exec_dir", "duration",
                      "direction", "traj_name", "vert_motion", "model_height",
                      "receptors", "system_type", "met_dir", "binary_path",
                      "folder_name", "config_list", "ascdata_list",
                      "recep_file_path", 'write_config_list',
                      'write_ascdata_list', 'get_receptor_values',
                      'get_daily_filenames', 'to_short_year', 'to_short_month',
                      'to_short_day', 'formatC', 'get_traj_output_filename',
                      'write_traj_control_file', 'to_null_dev',
                      'execute_on_system'),
                envir = environment())
  print('clusters initialized')
  traj.lst <- parLapply(cl, receptors, function(receptor) {
    
    inner_folder <- paste0('receptor', receptor)
    inner_dir <- file.path(exec_dir, inner_folder)
    if (!dir.exists(inner_dir)) {
      dir.create(inner_dir)
    }
    
    config_list %>% write_config_list(dir = file.path(inner_dir))
    ascdata_list %>% write_ascdata_list(dir = file.path(inner_dir))
    
    receptor_vals <-
      get_receptor_values(
        receptors_tbl = receptors_tbl,
        receptor_i = receptor
      )
    
    receptor_i <- receptor_vals$receptor
    lat_i <- receptor_vals$lat
    lon_i <- receptor_vals$lon
    height_i <- receptor_vals$height
    date_i <- receptor_vals$date
    hour_i <- receptor_vals$hour
    
    met_files <-
      get_daily_filenames(
        as.Date(date_i), duration, direction, suffix = "_hysplit.t00z.namsa"
      )
    
    #list_run_days <- days %>% as.character()
    
    start_year_GMT <- to_short_year(date_i)
    start_month_GMT <- to_short_month(date_i)
    start_day_GMT <- to_short_day(date_i)
    
    # Sort daily starting hours if given as
    # numeric values
    if (inherits(hour_i, "numeric")) {
      hour_i <- formatC(sort(hour_i), width = 2, flag = 0)
    }
    
    start_hour_GMT <- hour_i
    full_year_GMT <- as.character(year(as.Date(date_i)))
    
    output_filename <-
      get_traj_output_filename(
        traj_name = traj_name,
        site = receptor_i,
        direction = direction,
        year = start_year_GMT,
        month = start_month_GMT,
        day = start_day_GMT,
        hour = start_hour_GMT,
        lat = lat_i,
        lon = lon_i,
        height = height_i,
        duration = duration
      )
    
    
    # Write the CONTROL file
    write_traj_control_file(
      start_year_GMT = start_year_GMT,
      start_month_GMT = start_month_GMT,
      start_day_GMT = start_day_GMT,
      start_hour_GMT = start_hour_GMT,
      lat = lat_i,
      lon = lon_i,
      height = height_i,
      direction = direction,
      duration = duration,
      vert_motion = vert_motion,
      model_height = model_height,
      met_files = met_files,
      output_filename = output_filename,
      system_type = system_type,
      met_dir = met_dir,
      exec_dir = inner_dir
    )
    
    # The CONTROL file is now complete and in the
    # working directory, so, execute the model run
    sys_cmd <-
      paste0(
        "(cd \"",
        inner_dir,
        "\" && \"",
        binary_path,
        "\" ",
        to_null_dev(system_type = system_type),
        ")"
      )
    
    execute_on_system(sys_cmd, system_type = system_type)
    
    
    #recep_file_path <- file.path(exec_dir, folder_name)
    if (!dir.exists(recep_file_path)) {
      dir.create(path = recep_file_path, recursive = TRUE)
    }
    
    # Move files into the output folder
    trajectory_file <- file.path(inner_dir, output_filename)
    file.copy(
      from = trajectory_file,
      to = recep_file_path,
      copy.mode = TRUE
    )
    
    unlink(inner_dir, force = TRUE, recursive = TRUE)
  })
  stopCluster(cl)
  
  print('parallel computing complete')
  
  # For every set of coordinates, perform a set
  # of model runs
  # Obtain a trajectory data frame
  traj_tbl <-
    trajectory_read(output_folder = recep_file_path) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      #receptor = receptor_i,
      lat_i = lat_i,
      lon_i = lon_i,
      height_i = height_i
    )
  
  if (clean_up) {
    unlink(recep_file_path, recursive = TRUE, force = TRUE)
  }
  
  ensemble_tbl <-
    traj_tbl %>%
    dplyr::select(-c(year, month, day, hour)) %>%
    dplyr::select(
      receptor,
      hour_along,
      traj_dt,
      lat,
      lon,
      height,
      traj_dt_i,
      lat_i,
      lon_i,
      height_i,
      dplyr::everything()
    ) %>%
    dplyr::group_by(
      receptor, hour_along, traj_dt, traj_dt_i, lat_i, lon_i, height_i) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  if (direction == "forward") {
    
    ensemble_tbl <-
      ensemble_tbl %>%
      dplyr::arrange(receptor, traj_dt_i)
    
  } else {
    
    ensemble_tbl <-
      ensemble_tbl %>%
      dplyr::arrange(receptor, traj_dt_i, dplyr::desc(hour_along))
  }
  
  ensemble_tbl %>%
    dplyr::right_join(
      ensemble_tbl %>%
        dplyr::select(receptor, traj_dt_i, lat_i, lon_i, height_i) %>%
        dplyr::distinct() %>%
        dplyr::mutate(run = dplyr::row_number()),
      by = c("receptor", "traj_dt_i", "lat_i", "lon_i", "height_i")
    ) %>%
    dplyr::select(run, dplyr::everything())
  return(ensemble_tbl)
}


## Receptor should be the outpupt of as.list(get_receptor_values)
## Should be a list containing: receptor, lat, lon, height, date, hour

## Met filenames should be retrieved before function is run
# met_files <-
#   get_daily_filenames(
#     as.Date(date_i), duration, direction, suffix = "_hysplit.t00z.namsa"
#   )

