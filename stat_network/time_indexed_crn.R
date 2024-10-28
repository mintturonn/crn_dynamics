

# CRN matrices, three dimensions
crn_screen_pr <- array(runif(popsize * (params[["simlength"]]  + 1) * num_draws), 
                       dim = c(popsize, params[["simlength"]]  + 1, num_draws))

crn_pn_pr <- array(runif(popsize * (params[["simlength"]]  + 1) * num_draws), 
                   dim = c(popsize, params[["simlength"]]  + 1, num_draws))