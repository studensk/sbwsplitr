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
                               local_time = FALSE,
                               rdf_write = TRUE,
                               rdf_write_name = NULL) {
  
  
  config_list <-  list(KMSL = 0,
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
  met_file_list <- lapply(days, function(d) {
    get_daily_filenames(d, duration, direction, suffix = "_hysplit.t00z.namsa")
  })
  met_file_check <- unique(unlist(met_file_list))
  
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
  # if (is.null(traj_name)) {
  #   folder_name <- paste0("traj-", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
  # } else if (!is.null(traj_name)) {
  #   folder_name <- traj_name
  # }
  # recep_file_path <- file.path(exec_dir, folder_name)
  
  
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
      # dplyr::select(lat, lon) %>%
      # dplyr::distinct() %>%
      # dplyr::mutate(receptor = dplyr::row_number()) %>%
      # base::merge(run_df) %>%
      # dplyr::arrange(receptor) %>%
      dplyr::select(receptor, dplyr::everything())
  }
  
  else {
    receptors_tbl <- run_df %>%
      dplyr::as_tibble() %>%
      # dplyr::select(lat, lon) %>%
      # dplyr::distinct() %>%
      # dplyr::mutate(receptor = dplyr::row_number()) %>%
      # base::merge(run_df) %>%
      # dplyr::arrange(receptor, height, hour) %>%
      dplyr::select(receptor, dplyr::everything()) %>%
      dplyr::mutate(traj.name = paste('traj', date, receptor,
                                      height, hour, sep = '_'))
    
  }
  
  if (rdf_write) {
    if (is.null(rdf_write_name)) {rdf_write_name <- 'run_data'}
    write.csv(receptors_tbl, 
              paste0(file.path(exec_dir, rdf_write_name), '.csv'),
              row.names = FALSE)
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
  receptors <- unique(receptors_tbl$receptor)
  
  recep_file_path <- file.path(exec_dir, 'receptor_files')
  dir.create(recep_file_path)
  
  cores <- parallel::detectCores()
  max_clusters <- floor(cores*2/3)
  clusters <- pmin(max_clusters, max(receptors))
  print(paste0(clusters, ' clusters'))
  
  print('initialize clusters')
  cl <- makeCluster(clusters)
  # clusterEvalQ(cl, {
  #   library(tidyverse)
  # })
  clusterExport(cl, c("receptors_tbl", "exec_dir", "duration",
                      "direction", "traj_name", "vert_motion", "model_height",
                      "receptors", "system_type", "met_dir", "binary_path",
                      #"folder_name", 
                      'recep_file_path', 
                      "config_list", "ascdata_list",
                      "recep_file_path", 'write_config_list',
                      'write_ascdata_list', 'get_receptor_values',
                      'get_daily_filenames', 'to_short_year', 'to_short_month',
                      'to_short_day', 'formatC', 'get_traj_output_filename',
                      'write_traj_control_file', 'to_null_dev',
                      'execute_on_system', 'days', 'met_file_list',
                      'trajectory_read', 'tidy_grepl', 'tidy_gsub',
                      'execute_hysplit', 'get_os'),
                envir = environment())
  print('clusters initialized')
  
  for (day.ind in 1:length(days)) {
    clusterExport(cl, 'day.ind')
    
    clusterEvalQ(cl, {
      library(tidyverse)
      d <- days[day.ind]
      date_tbl <- subset(receptors_tbl, date == d)
      receptors <- unique(date_tbl$receptor)
      met_files <- met_file_list[[day.ind]]
    } )
    
    
    traj.lst <- parLapply(cl, receptors, function(rec) {
      
      inner_recep_file_path <- file.path(recep_file_path,
                                         paste0('receptor', rec))
      if (!dir.exists(inner_recep_file_path)) {
        dir.create(inner_recep_file_path)
      }
      
      rec_tbl <- subset(date_tbl, receptor == rec)
      
      inner_folder <- paste0('receptor', rec)
      inner_dir <- file.path(exec_dir, inner_folder)
      if (!dir.exists(inner_dir)) {
        dir.create(inner_dir)
      }
      
      config_list %>% write_config_list(dir = file.path(inner_dir))
      ascdata_list %>% write_ascdata_list(dir = file.path(inner_dir))
      
      for (tn in unique(rec_tbl$traj.name)) {
        receptor_vals <-
          get_receptor_values(
            receptors_tbl = rec_tbl,
            receptor_value = rec,
            traj_name = tn
          )
        
        execute_hysplit(receptor_vals = receptor_vals,
                        met_dir = met_dir,
                        met_files = met_files, 
                        directory = inner_dir,
                        inner_recep_file_path = inner_recep_file_path)
        
        
        # trajectory_file <- file.path(inner_dir, output_filename)
        # file.copy(
        #   from = trajectory_file,
        #   to = inner_recep_file_path,
        #   copy.mode = TRUE
        # )
      }
      
      unlink(inner_dir, force = TRUE, recursive = TRUE)
      
      if (day.ind == length(days)) {
        traj_tbl <-
          trajectory_read(output_folder = inner_recep_file_path) %>%
          dplyr::as_tibble() %>%
          dplyr::mutate(
            receptor = rec,
            lat_i = lat_i,
            lon_i = lon_i,
            height_i = height_i
          )

        write_file <- paste0('trajectories_receptor_', rec, '.csv')
        write_path <- file.path(recep_file_path, write_file)
        write.csv(traj_tbl, write_path, row.names = FALSE)
        
        unlink(inner_recep_file_path, recursive = TRUE)
      }
    
    })
    
  }
  # traj_tbl <-
  #   trajectory_read(output_folder = recep_file_path) %>%
  #   dplyr::as_tibble() %>%
  #   dplyr::mutate(
  #     #receptor = receptor_i,
  #     lat_i = lat_i,
  #     lon_i = lon_i,
  #     height_i = height_i
  #   )
  # unlink(recep_file_path, recursive = TRUE)
  stopCluster(cl)
  print('parallel computing complete')
  
  # For every set of coordinates, perform a set
  # of model runs
  # Obtain a trajectory data frame
  
  
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
  
  ensemble_tbl_final <- ensemble_tbl %>%
    dplyr::right_join(
      ensemble_tbl %>%
        dplyr::select(receptor, traj_dt_i, lat_i, lon_i, height_i) %>%
        dplyr::distinct() %>%
        dplyr::mutate(run = dplyr::row_number()),
      by = c("receptor", "traj_dt_i", "lat_i", "lon_i", "height_i")
    ) %>%
    dplyr::select(run, dplyr::everything()) 
  
  return(ensemble_tbl_final)
}

execute_hysplit <- function(receptor_vals = NULL,
                            traj_name = NULL,
                            receptor = 1,
                            date = NULL,
                            hour = NULL,
                            lat = NULL,
                            lon = NULL,
                            height = NULL,
                            duration = 9,
                            direction = 'forward',
                            model_height = 20000,
                            vert_motion = 0,
                            met_files = NULL,
                            met_dir = file.path(getwd(), 'meteorology'),
                            directory = NULL,
                            inner_recep_file_path = NULL) {
  
  if (!is.null(receptor_vals)) {
    receptor <- receptor_vals$receptor
    lat <- receptor_vals$lat
    lon <- receptor_vals$lon
    height <- receptor_vals$height
    date <- receptor_vals$date
    hour <- receptor_vals$hour
    traj_name <- receptor_vals$traj.name
  }
  
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
  
  if (substr(traj_name, 1, 4) != 'traj') {
    traj_name <- paste0('traj_', traj_name)
  }
  
  output_filename <-
    get_traj_output_filename(
      traj_name = traj_name,
      site = receptor,
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
    met_files = met_files,
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
      binary_path,
      "\" ",
      to_null_dev(system_type = system_type),
      ")"
    )
  
  execute_on_system(sys_cmd, system_type = system_type)
  
  trajectory_file <- file.path(directory, output_filename)
  file.copy(
    from = trajectory_file,
    to = inner_recep_file_path,
    copy.mode = TRUE
  )
  
}
