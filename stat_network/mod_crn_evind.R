
sims_pn_crn_ei <- function(init_inf, types, params, params_pn, tte_screen, tte_pn) {
  
  # Initialize matrices
  inf_state <- new_inf <- clearance <- diagn <- pn_index <- pn_partn <- matrix(0, length(init_inf), params[["simlength"]] + 1)
  inf_state[, 1] <- init_inf
  
  # List index cases' partner types
  partn <- apply(types, 1, function(x) which(x > 0))
  
  for (i in 2:(params[["simlength"]] + 1)) {
    for (m in 1:length(init_inf)) {
      
      # Infection acquisition
      if (inf_state[m, i - 1] == 0) {
        main_inf <- (rbinom(length(inf_state[partn[[m]], i - 1]), 1, params[["trnsm_pr"]])) == TRUE
        new_inf[m, i] <- ifelse(any(main_inf) > 0, 1, 0)
        inf_state[m, i] <- new_inf[m, i]
        
      } else {
        clearance[m, i] <- rbinom(1, 1, params[["clear_pr"]])
        
        # Diagnose at exact screening time
        if (length(tte_screen[[m]]) > 0 & i == tte_screen[[m]][1]) {
          diagn[m, i] <- 1
          tte_screen[[m]] <- tte_screen[[m]][-1]  # Remove used index
        }
        
        # Update infection state
        inf_state[m, i] <- ifelse(clearance[m, i] > 0 | diagn[m, i] > 0, 0, inf_state[m, i - 1])
      }
      
      # Partner notification as a binary event
      if (diagn[m, i] > 0 & params[["pn_pr"]] > 0) {
        if (tte_pn[[m]][i] == 1) {
          pn_index[m, i] <- 1
        }
        
        if (pn_index[m, i] > 0) {
          # Treat all partners
          inf_state[partn[[m]], i] <- 0
          pn_partn[partn[[m]], i] <- 3
          diagn[partn[[m]], i] <- ifelse(inf_state[partn[[m]], i - 1] > 0, 3, 0)
        }
      }
    }
  }
  
  return(list(inf = inf_state, clear = clearance, diagn = diagn, pn_partn = pn_partn, new_inf = new_inf))
}