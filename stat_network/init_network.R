
####################################
# initialize the network

id <- paste("M", 1:popsize, sep="")
linkm <- matrix(0, popsize, popsize, dimnames = list(id, id))

# everyone has at least one tie and most have 1-5 ties
for (i in 1:popsize) {
  
  num_connect <- sample(1:8, 1, prob = c(0.4, 0.44, 0.1, 0.01, 0.01, 0.01, 0.01, 0.02))
  partners <- sample(setdiff(1:popsize, i), num_connect, replace = FALSE)
  linkm[i, partners] <- 1
}

# symmetric
linkm <- linkm + t(linkm)
linkm[linkm > 1] <- 1

diag(linkm) <- 0  

# test symmetry 
linkm[lower.tri(linkm)] <- t(linkm)[lower.tri(linkm)]
