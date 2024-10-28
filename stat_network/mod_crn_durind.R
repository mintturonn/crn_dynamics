

sims_pn_crn_di <- function(init_inf, types,  params, params_pn, tte_screen, tte_pn) {
  
  inf_state <- new_inf <- clearance <- diagn <- pn_index <- inf_dur <-  pn_partn <- matrix(0, length(init_inf), params[["simlength"]] +1)
  inf_state[,1] <- init_inf
  pn_coverage <-  vector()
  
  # list index cases' partner types
  partn <- apply(types, 1, function (x) which(x>0))
  
  
  for (i in 2:(params[["simlength"]] +1)){
    for (m in 1:length(init_inf) ){
      
      # infection acquisition
      if ( inf_state[m,i-1]==0 ) {
        
        # as long as everyone has >0 partners, this works
        main_inf <- (rbinom(length(inf_state[partn[[m]],i-1]) ,1, params[["trnsm_pr"]]))==TRUE 
        
        new_inf[m,i] <- ifelse(any(main_inf) >0, 1, 0)
        
        inf_state[m,i] <- new_inf[m,i]
        # infection clearance  
      }else{
        inf_dur[m, i] <- inf_dur[m, i - 1] + 1
        clearance[m,i] <- rbinom(1,1, params[["clear_pr"]])

        # Diagnose infection after X time steps
        diagn[m, i] <- ifelse(inf_dur[m, i] == tte_screen[[m]][1], 1, 0)
        
        if (clearance[m, i] > 0 | diagn[m, i] > 0) {
            inf_state[m, i] <- 0
            inf_dur[m, i] <- 0
            tte_screen[[m]] <- tte_screen[[m]][-1]
        } else {
          inf_state[m, i] <- 1
        }
      }
        
      # partner notification
      if ( diagn[m,i] >0 &  params[["pn_pr"]]>0  ) {
        
        # intervention increase
        # ifelse(i>209, params[["pn_pr"]] <- params_pn, params[["pn_pr"]])
        
        pn_index[m,i] <- tte_pn[[m]][[1]]
        tte_pn[[m]] <- tte_pn[[m]][-1]
        
        if (  pn_index[m,i] > 0) {
          
          # all partners treated
          inf_state[partn[[m]],i] <- 0
          pn_partn[partn[[m]],i] <- 3
          diagn[partn[[m]],i] <- ifelse(inf_state[partn[[m]],i-1]>0, 3, 0)
          
        }
      }
    }
  }
  
  return(list(inf = inf_state, clear = clearance, diagn = diagn, pn_partn=pn_partn, new_inf=new_inf))
}



