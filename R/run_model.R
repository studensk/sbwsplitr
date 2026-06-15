#' Run the model
#'
#' Run either a trajectory model or a dispersion model, depending on the class
#' of the splitr modeling object.
#' 
#' @param model A splitr modeling object.
#' @export
run_model <- function(model, particle_max = '100000') {
  
  if (inherits(model, "trajectory_model")) {
    
    traj_df <- 
      hysplit_trajectory(
        lat = model$lat,
        lon = model$lon,
        height = model$height,
        duration = ifelse(is.null(model$duration), 24, model$duration),
        days = model$days,
        daily_hours = model$daily_hours,
        direction = ifelse(is.null(model$direction),
                           "forward", model$direction),
        met_type = ifelse(is.null(model$met_type),
                          "reanalysis", model$met_type),
        vert_motion = ifelse(is.null(model$vert_motion),
                             0, model$vert_motion),
        model_height = ifelse(is.null(model$model_height),
                              20000, model$model_height),
        extended_met = TRUE,
        traj_name = model$traj_name,
        exec_dir = model$exec_dir,
        met_dir = model$met_dir,
        binary_path = model$binary_path,
        softrun = model$softrun,
        clean_up = model$clean_up
      )
    
    model$traj_df <- traj_df
    
    return(model)
  }
  
  if (inherits(model, "dispersion_model")) {
    
    # Get time window for observations
    start_day <- model$start_time %>% lubridate::floor_date()
    start_hour <- model$start_time %>% lubridate::hour() 
    duration <- as.numeric(difftime(model$end_time, model$start_time, units = "hours"))
    
    # Get ith source parameters
    # lat <- model$geo_df[["lat"]]
    # lon <- model$geo_df[["lon"]]
    # height <- model$geo_df[i, ][["height"]]
    
    # release_start <- model$sources[["release_start"]]
    # release_end <- model$sources[["release_end"]]
    # 
    # rate <- model$sources[["rate"]]
    
    species_list <- 
      model$sources %>% 
      #dplyr::select(-c(lat, lon, height)) %>%
      as.list()
    
    nsources <- nrow(model$geo_df)
    
    disp_df <- 
      hysplit_dispersion(
        # lat = lat,
        # lon = lon,
        # height = height,
        geo_df = model$geo_df,
        start_day = start_day,
        start_hour = start_hour,
        duration = duration,
        direction = "forward",
        met_type = model$met_type,
        vert_motion = model$vert_motion,
        model_height = model$model_height,
        # particle_num = 10*nsources, #2500,
        # particle_max = 100*nsources, #10000,
        particle_num = ifelse(is.null(model$particle_num), '-10', 
                              model$particle_num),
        particle_max = particle_max,
        species = species_list,
        exec_dir = model$exec_dir,
        met_dir = model$met_dir,
        binary_path = model$binary_path,
        binary_name = model$binary_name,
        softrun = model$softrun,
        clean_up = model$clean_up
      )
    
    model$disp_df <- disp_df
    return(model) 
  }
}
