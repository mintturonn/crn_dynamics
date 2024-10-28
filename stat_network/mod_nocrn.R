

sims_pn_nocrn <- function(init_inf, types,  params, params_pn) {
  
  inf_state <- new_inf <- clearance <- diagn <- pn_index <-  pn_partn <-  matrix(0, length(init_inf), params[["simlength"]] +1)
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
          clearance[m,i] <- rbinom(1,1, params[["clear_pr"]])
          
          # someone can be diagnosed as an index case and receive PN in the same time step 
          diagn[m,i]     <- ifelse(rbinom(1,1, params[["screen_pr"]])==1, 1, 0)  
          
          inf_state[m,i] <- ifelse(clearance[m,i]>0 | diagn[m,i]>0 , 0, inf_state[m,i-1])
        }
        
        # partner notification
        if ((diagn[m,i] >0 &  params[["pn_pr"]]>0 )>0 ) {
          
          # intervention increase
         # ifelse(i>209, params[["pn_pr"]] <- params_pn, params[["pn_pr"]])
          
          pn_index[m,i] <- ifelse(rbinom(1,1, params[["pn_pr"]])==1, 1, 0)  
          
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



