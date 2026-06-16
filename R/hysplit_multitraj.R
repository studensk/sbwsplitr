#' Conduct extract temperature data from multiple starting points in HYSPLIT
#' 
#' 
#' @inheritParams hysplit_trajectory
#'   
#' @export

hysplit_multitraj <- function(full_df,
                              direction = "forward",
                              met_type = "nam12",
                              vert_motion = 0,
                              model_height = 20000,
                              extended_met = TRUE,
                              temp_only = TRUE,
                              traj_name = NULL,
                              met_dir = NULL,
                              exec_dir = NULL,
                              system_type = 'win') {
  
  unlink(list.files(file.path(exec_dir), full.names = TRUE))
  
  disp.dh <- full_df |>
    select(date, hour) |>
    unique()
  
  download_met_files(met_type, unique(disp.dh$date), 
                     duration = 0, 
                     direction = direction, 
                     met_dir = met_dir)
  
  binary_path <- 
    set_binary_path(
      binary_path = NULL,
      binary_name = "hyts_std"
    )
  
  for (d in seq_len(nrow(disp.dh))) {
    
    hr <- disp.dh$hour[d]
    dat <- disp.dh$date[d]
    
    sub.df <- full_df |>
      subset(date == dat & hour == hr) |>
      select(lat, lon, height)
    
    hysplit_trajectory(
      geo_df = sub.df,
      duration = 0,
      date = dat,
      hour = hr,
      clean_up = FALSE,
      exec_dir = exec_dir,
      met_dir = met_dir,
      met_type = met_type,
      extended_met = extended_met,
      temp_only = temp_only
    )
    
    sys_cmd <- 
      paste0(
        "(cd \"",
        exec_dir,
        "\" && \"",
        binary_path,
        "\" ",
        to_null_dev(system_type = system_type),
        ")"
      )
    shell(sys_cmd)
  }
  
  
  standard_col_names <- 
    c(
      "year", "month", "day", "hour", 
      "lat", "lon", "height", "temperature"
    )
  
  traj.files <- list.files('trajectory_test', pattern = '^traj--')
  line.lst <- lapply(traj.files, function(file) {
    file.name <- paste0('trajectory_test/', file)
    file_lines <- readLines(file.name, encoding = "UTF-8", skipNul = TRUE)
    file_one_line <- readr::read_file(file.name)
    header_line <-
      file_lines %>%
      vapply(
        FUN.VALUE = logical(1),
        USE.NAMES = FALSE,
        function(x) tidy_grepl(x, "PRESSURE")
      ) %>%
      which()
    file_lines
    file_one_line
    header_line
    file_lines_data <-
      file_lines[(header_line + 1):(length(file_lines))] %>%
      tidy_gsub("\\s\\s*", " ") %>%
      tidy_gsub("^ ", "")
    return(file_lines_data)
  })
  all.lines <- unlist(line.lst)
  
  traj_tbl <- 
    all.lines %>%
    strsplit("\\s+") %>%
    lapply(
      FUN = function(x) {
        x[c(3:6, 10:12, 14)] %>%
          as.numeric() %>%
          stats::setNames(standard_col_names) %>%
          as.list() %>%
          dplyr::as_tibble()
      }
    ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(temperature = temperature - 273.15) %>%
    merge(full_df) 
  
  return(traj_tbl)
}


