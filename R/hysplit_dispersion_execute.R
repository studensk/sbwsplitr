#' Conduct HYSPLIT dispersion runs
#' 
#' The function executes single/multiple forward or backward HYSPLIT dispersion
#' runs using specified meteorological datasets.
#' 
#' @param disp_setup List output from hysplit_dispersion.R
#' @export

hysplit_dispersion_execute <- function(disp_setup) {
  # The CONTROL file is now complete and in the
  # working directory, so, execute the model run
  output <- with(disp_setup, {
    sys_cmd <- 
      paste0(
        "(cd \"",
        exec_dir,
        "\" && \"",
        hycs_std_binary_path,
        "\" ",
        to_null_dev(system_type = system_type),
        ")"
      )
    
    execute_on_system(sys_cmd, system_type = system_type)
    
    # Extract the particle positions at every hour
    sys_cmd <- 
      paste0(
        "(cd \"",
        exec_dir,
        "\" && \"",
        par2asc_binary_path,
        "\" -iPARDUMP -opardump_output.txt",
        to_null_dev(system_type = system_type),
        ")"
      )
    
    execute_on_system(sys_cmd, system_type = system_type)
    dispersion_tbl <- parse_pardump(file.path(exec_dir, 'pardump_output.txt'))
    
    if (clean_up) {
      unlink(file.path(exec_dir, list.files(path = exec_dir,
                                            pattern = "^.*$")), force = TRUE)
    }
    
    dispersion_tbl
  })
  return(output)
}