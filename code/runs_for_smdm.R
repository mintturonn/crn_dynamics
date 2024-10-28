

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
source(here('stat_network/mod_crn_evind.R'))
source(here('stat_network/mod_crn_timeind.R'))

#################


popsize <- 1000
# 
set.seed(4321)
source(here('stat_network/init_network.R'))

hist(rowSums(linkm), breaks=100, ylab="Frequnecy", xlab="Number of partners", main="")
quantile(rowSums(linkm))

params <- NULL
params[["timehorizon"]] <- 2
params[["timestep"]] <- 52
params[["simlength"]] <- params[["timehorizon"]]*params[["timestep"]]
params[["popsize"]] <- popsize
params[["clear_pr"]] <- 0
params[["trnsm_pr"]] <- 0.01
params[["screen_pr"]] <- 0.08
params[["pn_pr"]] <-  0#0.5

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

ort_p <- ort_i <- numeric()

t0 <- Sys.time() 

for (i in 1:num_draws){
  
  ort <-  sims_pn_crn_ti(init_inf, linkm, params, params_pn, crn_screen_pr[,,i], crn_pn_pr[,,i])
  ort_p <-   rbind(ort_p, colSums(ort$inf)) 
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
  theme_minimal() + ylim(c(0, NA)) -> pt

ggsave(here('output/timeind_prev_nopn.png'), plot = pt, width = 4, height = 4, dpi = 300)


#################

source(here('stat_network/stat_network/duration_indexed_crn.R'))

# duration-indexed CRN

set.seed(321)
# with 100 draws, this takes a while 
tte <- generate_tte(popsize=popsize, num_draws=num_draws, params=params)


ord_p <- ord_i <- numeric()

t0 <- Sys.time() 

for (i in 1:num_draws){
  
  ord <-  sims_pn_crn_di(init_inf, linkm, params, params_pn, tte$tte_screen[[i]], tte$tte_pn[[i]])
  ord_p <-   rbind(ord_p, colSums(ord$inf)) 
  ord_i <-   rbind(ord_i, colSums(ord$new_inf)) 
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
  theme_minimal() + ylim(c(0, NA)) -> pd

ggsave(here('output/durind_prev_nopn.png'), plot = pd, width = 4, height = 4, dpi = 300)


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

ggsave(here('output/eventind_prev_nopn.png'), plot = pe, width = 4, height = 4, dpi = 300)



#################

set.seed(2)

or_p <- or_i <- numeric()

t0 <- Sys.time() 

for (i in 1:100){
  
  or <-  sims_pn_nocrn(init_inf, linkm, params, params_pn)
  or_p <-   rbind(or_p, colSums(or$inf)) 
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
  ylim(c(0, NA)) -> ps
  
ggsave(here('output/stochastic_prev_nopn.png'), plot = ps, width = 4, height = 4, dpi = 300)



