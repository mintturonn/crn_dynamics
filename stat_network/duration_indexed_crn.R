
generate_tte <- function(popsize, num_draws, params) {
  tte_screen <- replicate(num_draws, {
    lapply(1:popsize, function(x) {
      steps <- numeric()
      for (j in 1:100) {
        n_steps <- 0
        detected <- 0
        while (detected == 0) {
          n_steps <- n_steps + 1
          detected <- rbinom(1, 1, params[["screen_pr"]])
        }
        steps <- c(steps, n_steps)
      }
      steps
    })
  }, simplify = FALSE)
  
  tte_pn <- replicate(num_draws, {
    lapply(1:popsize, function(x) { notified <- rbinom(100, 1, params[["pn_pr"]]) })
  }, simplify = FALSE)
  
  tte_nopn <- replicate(num_draws, {
    lapply(1:popsize, function(x) { notified <- rbinom(100, 1, 0) })
  }, simplify = FALSE)
  
  tte_allpn <- replicate(num_draws, {
    lapply(1:popsize, function(x) { notified <- rbinom(100, 1, 1) })
  }, simplify = FALSE)
  
  list(tte_screen = tte_screen, tte_pn = tte_pn, tte_nopn = tte_nopn, tte_allpn=tte_allpn)
}