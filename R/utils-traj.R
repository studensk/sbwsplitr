#' @noRd
write_traj_control_file <- function(start_year_GMT,
                                    start_month_GMT,
                                    start_day_GMT,
                                    start_hour_GMT,
                                    # lat,
                                    # lon,
                                    # height,
                                    geo_df,
                                    direction,
                                    duration,
                                    vert_motion,
                                    model_height,
                                    met_files,
                                    output_filename,
                                    system_type,
                                    met_dir,
                                    exec_dir) {
  
  npts <- nrow(geo_df)
  geo_string <- geo_df |>
    mutate(string = paste0(lat, " ", lon, " ", height, "\n")) |>
    select(string) |>
    unlist() |>
    paste0(collapse = '')
  
  paste0(
    start_year_GMT, " ", start_month_GMT, " ",
    start_day_GMT, " ", start_hour_GMT, "\n",
    npts, "\n",
    #lat, " ", lon, " ", height, "\n",
    geo_string,
    ifelse(direction == "backward", "-", ""), duration, "\n",
    vert_motion, "\n",
    model_height, "\n",
    length(met_files), "\n",
    paste0(met_dir, "/\n", met_files, collapse = "\n"), "\n",
    exec_dir, "/\n",
    output_filename, "\n"
  ) %>%
    cat(file = file.path(exec_dir, "CONTROL"), sep = "", append = FALSE)
}
