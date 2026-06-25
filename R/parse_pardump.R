parse_pardump <- function(file) {
  
  con <- file(file, "r")
  on.exit(close(con))
  
  particle_id <- numeric()
  page_hr     <- numeric()
  lat         <- numeric()
  lon         <- numeric()
  height      <- numeric()
  
  n <- 0
  
  repeat {
    
    line <- readLines(con, n = 1)
    
    if (length(line) == 0)
      break
    
    line <- trimws(line)
    
    # Skip description blocks
    if (startsWith(line, "Header Record:")) {
      readLines(con, n = 3)
      line <- trimws(readLines(con, n = 1))
    }
    
    hdr <- scan(text = line, quiet = TRUE)
    numpar <- hdr[1]
    
    for (p in seq_len(numpar)) {
      
      # skip mass
      readLines(con, n = 1)
      
      pos <- scan(text = readLines(con, n = 1), quiet = TRUE)
      
      meta <- scan(text = readLines(con, n = 1), quiet = TRUE)
      
      n <- n + 1
      
      particle_id[n] <- meta[5]
      page_hr[n]     <- meta[1] / 60
      lat[n]         <- pos[1]
      lon[n]         <- pos[2]
      height[n]      <- pos[3]
    }
  }
  
  data.frame(
    particle_id = particle_id,
    hour_along     = page_hr,
    lat         = lat,
    lon         = lon,
    height      = height
  )
}
