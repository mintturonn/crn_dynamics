

generate_tte2 <- function(popsize, num_draws, params) {
  
  generate_screen_indices <- function(pr, simlength) {
    steps <- numeric()
    window_starts <- seq(1, simlength, length.out = 10)
    
    for (start in window_starts) {
      event_time <- which(runif(10) < pr)[1]
      if (!is.na(event_time)) {
        steps <- c(steps, round(start + event_time - 1))
      }
    }
    
    steps
  }
  
  # Retain binary logic for partner notification
  generate_pn_events <- function(pr, simlength) {
    rbinom(simlength, 1, pr)
  }
  
  tte_screen <- replicate(num_draws, {
    lapply(1:popsize, function(x) {
      generate_screen_indices(params[["screen_pr"]], params[["simlength"]])
    })
  }, simplify = FALSE)
  
  tte_pn <- replicate(num_draws, {
    lapply(1:popsize, function(x) {
      generate_pn_events(params[["pn_pr"]], params[["simlength"]])
    })
  }, simplify = FALSE)
  
  list(tte_screen = tte_screen, tte_pn = tte_pn)
}