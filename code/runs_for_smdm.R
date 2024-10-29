

library(here)
library(dplyr)
library(reshape2)
library(tidyverse)
library(readxl)


# rm(list = ls())
# .rs.restartR()
###################################

# transmission models 
source(here('stat_network/mod_nocrn.R'))
source(here('stat_network/mod_crn_timeind.R'))
source(here('stat_network/mod_crn_durind.R'))
source(here('stat_network/mod_crn_durind_fixed.R'))

# this does not work: 
source(here('stat_network/mod_crn_evind.R'))

#################


popsize <- 1000
# 
set.seed(4321)
source(here('stat_network/init_network.R'))

hist(rowSums(linkm), breaks=100, ylab="Frequnecy", xlab="Number of contacts", main="")
quantile(rowSums(linkm))

params <- NULL
params[["timehorizon"]] <- 4
params[["timestep"]] <- 52
params[["simlength"]] <- params[["timehorizon"]]*params[["timestep"]]
params[["popsize"]] <- popsize
params[["clear_pr"]] <- 0
params[["trnsm_pr"]] <- 0.06
params[["screen_pr"]] <- 0.1
params[["pn_pr"]] <-  0.5

params_pn <- params[["pn_pr"]] 

#################

# initial infection 

set.seed(1)
init_inf <- rbinom(popsize, 1, 0.2) 


#################

# time-indexed CRN

set.seed(123)  
num_draws=100
source(here('stat_network/time_indexed_crn.R'))

# i=1

ort_p <-  ort_ip <- ort_i <- numeric()

t0 <- Sys.time() 

for (i in 1:num_draws){
  
  ort <-  sims_pn_crn_ti(init_inf, linkm, params, params_pn, crn_screen_pr[,,i], crn_pn_pr[,,i])
  ort_p <-   rbind(ort_p, colSums(ort$inf)) 
  ort_ip <-   rbind(ort_ip, rowSums(ort$inf)) 
  ort_i <-   rbind(ort_i, colSums(ort$new_inf)) 
}

Sys.time()-t0


ort_i %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", time))) %>%
  ggplot( aes(x = Time, y = Value/popsize, group = sim)) +
  geom_line(alpha=0.2) +
  labs(title = "incidence",
       x = "Time Step",
       y = "Value") +
  theme_minimal() + ylim(c(0, NA))

ort_p %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "Time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", Time))) %>%
  group_by(Time) %>%
  mutate(medianvalue = median(Value)) %>%  # Calculate the median for each time step
  ggplot(aes(x = Time, y = 100 * Value / popsize, group = sim)) +
  geom_line(alpha = 0.1) +
  geom_line(aes(y = 100 * medianvalue / popsize), color = "red", size = 1) +  
  labs(title = "time-index: prevalence",
       x = "Time Step",
       y = "%") +
  theme_minimal() +  ylim(c(0, 60)) -> pt

ggsave(here('output/timeind_prev_pn.png'), plot = pt, width = 4, height = 4, dpi = 300)

quantile(ort_ip, probs = c(0.5, 0.025, 0.975))/52


#################

source(here('stat_network/duration_indexed_crn.R')) 

# duration-indexed CRN

# set.seed(123)
# with 100 draws, this takes a while 
tte <- generate_tte(popsize=popsize, num_draws=num_draws, params=params)


ord_p <- ord_ip <- ord_i <- ord_d <- ord_pn <- numeric()

t0 <- Sys.time() 

for (i in 1:num_draws){
  
  ord <-  sims_pn_crn_di(init_inf, linkm, params, params_pn, tte$tte_screen[[i]], tte$tte_pn[[i]])
  ord_p <-   rbind(ord_p, colSums(ord$inf)) 
  ord_ip <-   rbind(ord_ip, rowSums(ord$inf)) 
  ord_i <-   rbind(ord_i, colSums(ord$new_inf)) 
  ord_d <-   rbind(ord_i, colSums(ord$diagn==1)) 
  ord_pn <-   rbind(ord_i, colSums(ord$diagn==3)) 
}

Sys.time()-t0

ord_i %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", time))) %>%
  ggplot( aes(x = Time, y = Value/popsize, group = sim)) +
  geom_line(alpha=0.2) +
  labs(title = "incidence",
       x = "Time Step",
       y = "Value") +
  theme_minimal() + ylim(c(0, NA))

ord_p %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "Time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", Time))) %>%
  group_by(Time) %>%
  mutate(medianvalue = median(Value)) %>%  # Calculate the median for each time step
  ggplot(aes(x = Time, y = 100 * Value / popsize, group = sim)) +
  geom_line(alpha = 0.1) +
  geom_line(aes(y = 100 * medianvalue / popsize), color = "red", size = 1) +  
  labs(title = "duration-index: prevalence",
       x = "Time Step",
       y = "%") +
  theme_minimal() + ylim(c(0, 60))  -> pd

ggsave(here('output/durind_prev_pn.png'), plot = pd, width = 4, height = 4, dpi = 300)

quantile(ord_ip, probs = c(0.5, 0.025, 0.975))/52

############################


ord2_p <- ord2_i <- ord2_d <- ord2_pn <- ord2_ip <- numeric()

t0 <- Sys.time() 

for (i in 1:num_draws){
  
  ord2 <-  sims_pn_crn_di_fix(init_inf, linkm, params, params_pn, tte$tte_screen[[i]], tte$tte_pn[[i]])
  ord2_p <-   rbind(ord2_p, colSums(ord2$inf)) 
  ord2_ip <-   rbind(ord2_ip, rowSums(ord2$inf)) 
  ord2_i <-   rbind(ord2_i, colSums(ord2$new_inf)) 
  ord2_d <-   rbind(ord2_i, colSums(ord2$diagn==1)) 
  ord2_pn <-   rbind(ord2_i, colSums(ord2$diagn==3)) 
}

Sys.time()-t0

ord_i %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", time))) %>%
  ggplot( aes(x = Time, y = Value/popsize, group = sim)) +
  geom_line(alpha=0.2) +
  labs(title = "incidence",
       x = "Time Step",
       y = "Value") +
  theme_minimal() + ylim(c(0, NA))

ord2_p %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "Time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", Time))) %>%
  group_by(Time) %>%
  mutate(medianvalue = median(Value)) %>%  # Calculate the median for each time step
  ggplot(aes(x = Time, y = 100 * Value / popsize, group = sim)) +
  geom_line(alpha = 0.1) +
  geom_line(aes(y = 100 * medianvalue / popsize), color = "red", size = 1) +  
  labs(title = "duration-index: prevalence",
       x = "Time Step",
       y = "%") +
  theme_minimal() + ylim(c(0, 60))  -> pd

ggsave(here('output/durind_fixed_prev_pn.png'), plot = pd, width = 4, height = 4, dpi = 300)

# this is taking the crude summary statistic (all people represented num_draws times in the model)
quantile(ord2_ip, probs = c(0.5, 0.025, 0.975))/52


#############################

set.seed(2)

or_p <- or_i <- or_ip <- numeric()

t0 <- Sys.time() 

for (i in 1:200){
  
  or <-  sims_pn_nocrn(init_inf, linkm, params, params_pn)
  or_p <-   rbind(or_p, colSums(or$inf)) 
  or_ip <-   rbind(or_ip, rowSums(or$inf)) 
  or_i <-   rbind(or_i, colSums(or$new_inf)) 
}

Sys.time()-t0

plot(0:params[["simlength"]], colSums(or$inf))

or_i %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", time))) %>%
  ggplot( aes(x = Time, y = Value/popsize, group = sim)) +
  geom_line(alpha=0.2) +
  labs(title = "incidence",
       x = "Time Step",
       y = "Value") +
  theme_minimal() + ylim(c(0, NA))

or_p %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "Time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", Time))) %>%
  group_by(Time) %>%
  mutate(medianvalue = median(Value)) %>%  # Calculate the median for each time step
  ggplot(aes(x = Time, y = 100 * Value / popsize, group = sim)) +
  geom_line(alpha = 0.1) +
  geom_line(aes(y = 100 * medianvalue / popsize), color = "red", size = 1) +  
  labs(title = "Stochastic, prevalence",
       x = "Time Step",
       y = "%") +
  theme_minimal() +
  ylim(c(0, 60)) -> ps
  
ggsave(here('output/stochastic_prev_nopn.png'), plot = ps, width = 4, height = 4, dpi = 300)

quantile(or_ip, probs = c(0.5, 0.025, 0.975))/52


#################
# event-indexed CRN

source(here('stat_network/stat_network/event_indexed_crn.R'))

set.seed(21)
# with 100 draws, this takes a while 
tte2 <- generate_tte2(popsize=popsize, num_draws=num_draws, params=params)


ore_p <- ore_i <- numeric()

t0 <- Sys.time() 

for (i in 1:num_draws){
  
  ore <-  sims_pn_crn_ei(init_inf, linkm, params, params_pn, tte2$tte_screen[[i]], tte2$tte_pn[[i]])
  ore_p <-   rbind(ore_p, colSums(ore$inf)) 
  ore_i <-   rbind(ore_i, colSums(ore$new_inf)) 
}

Sys.time()-t0

ore_i %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "time", values_to = "Value") %>%
  mutate(Time = as.numeric(sub("V", "", time))) %>%
  ggplot( aes(x = Time, y = Value/popsize, group = sim)) +
  geom_line(alpha=0.2) +
  labs(title = "incidence",
       x = "Time Step",
       y = "Value") +
  theme_minimal() + ylim(c(0, NA))

ore_p %>%
  as_tibble() %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "Time", values_to = "Value") %>%
  group_by(Time) %>%
  mutate(medianvalue = median(Value)) %>%  # Calculate the median for each time step
  ggplot(aes(x = Time, y = 100 * Value / popsize, group = sim)) +
  geom_line(alpha = 0.1) +
  geom_line(aes(y = 100 * medianvalue / popsize), color = "red", size = 1) +  
  labs(title = "event-index: Prevalence",
       x = "Time Step",
       y = "%") +
  theme_minimal() + ylim(c(0, NA)) -> pe

ggsave(here('output/eventind_prev_pn0.5.png'), plot = pe, width = 4, height = 4, dpi = 300)



#################

