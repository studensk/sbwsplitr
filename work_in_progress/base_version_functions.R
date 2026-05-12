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

plot_start_points <- function(run_df,
                              pt_size = 1.5) {
  
  if (all(names(run_df) != 'date')) {
    run_df$date <- run_df$date_i
  }
  ll.distinct <- run_df %>%
    dplyr::mutate(date = as.Date(date),
                  year = year(date),
                  date_i = date) %>%
    #dplyr::select(lon, lat, year) %>%
    dplyr::select(lon, lat, date_i) %>%
    dplyr::distinct()
  g <- ggplot() +
    geom_point(data = ll.distinct,  aes(x = lon, y = lat), size = pt_size) + 
    #facet_wrap(vars(year)) +
    labs(x = 'Longitude', y = 'Latitude') +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 12, 
                                      margin = margin(t = 10)),
          axis.title.y = element_text(size = 12, 
                                      margin = margin(r = 10)))
  return(g)
  
}
plot_trajectories <- function(traj_df,
                              line_transparency = 1,
                              line_width = 1.25, 
                              origin_pt_size = 2,
                              traj_pt_size = 1.5, 
                              title = NULL,
                              facet_dates = FALSE) {
  
  traj_df <- traj_df %>%
    mutate(date_i = as.Date(date_i))
  start_pts <- subset(traj_df, hour_along == 0)
  g <- plot_start_points(start_pts,
                         pt_size = origin_pt_size)
  
  g.traj <- g + 
    geom_path(data = traj_df, aes(x = lon, y = lat,
                                  group = traj_id_full, col = factor(height_i)),
              alpha = line_transparency, linewidth = line_width) +
    geom_point(data = traj_df, aes(x = lon, y = lat,
                                   col = factor(height_i)), size = traj_pt_size) +
    geom_point(data = start_pts, aes(x = lon, y = lat), size = traj_pt_size) +
    labs(col = 'Starting \nHeight (m)', title = title) +
    theme(legend.title = element_text(hjust = 0.5))
  if (facet_dates) {
    g.traj <- g.traj + facet_wrap(vars(date_i))
  }
  return(g.traj)
}

multiple_trajectories <- function(run_df = NULL,
                                  date = NULL,
                                  hour = NULL,
                                  lat = NULL,
                                  lon = NULL,
                                  height = NULL,
                                  duration = 9,
                                  direction = 'forward',
                                  model_height = 20000,
                                  vert_motion = 0,
                                  #met_dir = file.path(getwd(), 'meteorology'),
                                  met_dir = getwd(),
                                  exec_dir = getwd(), 
                                  binary_path = NULL,
                                  csv_folder = 'traj_output',
                                  clean_up = TRUE,
                                  rdf_write = TRUE,
                                  rdf_write_name = 'run_data', 
                                  traj_write = TRUE,
                                  traj_write_name = 'trajectories_final',
                                  plot = FALSE, 
                                  vbug = 2.5, ...) {
  
  if (is.null(run_df)) {
    run_df <- make_run_df(lat, lon, height, date, hour)
  }
  
  
  # If the execution dir isn't specified, use the working directory
  if (!dir.exists(exec_dir)) {dir.create(exec_dir)}
  
  output_path <- file.path(exec_dir, csv_folder)
  if (!dir.exists(output_path)) {dir.create(output_path)}
  
  # If the meteorology dir isn't specified, use the working directory
  # if (!dir.exists(met_dir)) {dir.create(met_dir)}
  
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
  
  if (rdf_write) {
    if (is.null(rdf_write_name)) {rdf_write_name <- 'run_data'}
    write.csv(run_tbl, 
              paste0(file.path(exec_dir, rdf_write_name), '.csv'),
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
  if (traj_write) {
    # write.csv(traj.df, file.path(exec_dir, paste0(traj_write_name, '.csv')), 
    #           row.names = FALSE)
    write.csv(run.df.mg, file.path(exec_dir, paste0(traj_write_name, '.csv')), 
              row.names = FALSE)
  }
  
  if (clean_up) {unlink(output_path, recursive = TRUE, force = TRUE)}
  
  if (plot) {
    print(plot_trajectories(run.df.mg, ...))
  }
  return(run.df.mg)
}

make_run_df <- function(lat, lon, height, date, hour) {
  ll.df <- data.frame(lat, lon)
  df <- expand.grid(height, hour, date)
  names(df) <- c('height', 'hour', 'date')
  mg.df <- merge(df, ll.df)
  return(mg.df)
}
