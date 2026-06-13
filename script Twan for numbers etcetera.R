set.seed(42)

## Load packages
library(odin2)
library(dust2)
library(monty)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)

root<-here()

# Call model script5
source(file.path(root,"model with screening delay compartment.R"))


# Create model object
sys <- dust_system_create(uli_ode, 
                          pars = list(beta=3,beta_p=22))

# Set default initial state of the model 
dust_system_set_state_initial(sys)

# Run simulation from time zero
t <- seq(0, 1000)
y <- dust_system_simulate(sys, t)

# Get model state at the latest point for future simulations
init_state0<-dust_system_state(sys)

y <- dust_unpack_state(sys, y)

names(init_state0)<-paste0(names(y))

totalpop <- y$Is + y$Ia + y$U + y$L + y$R + y$Is_p + y$Ia_p + y$U_p + y$L_p + y$R_p + y$Tr + y$Tr_p + y$Sd_p + y$Sd

prison_pop<-y$Is_p + y$Ia_p + y$U_p + y$L_p + y$R_p + y$Tr_p + y$Sd_p

inc<-c(0,diff(y$incidence)*1e5)

inc_p<-c(0,diff(y$incidence_p))/prison_pop*1e5

prev_p<-1e2*(y$Is_p + y$Ia_p + y$Tr_p)/prison_pop 


# Plot baseline trajectories
plot(t ,inc, type = "l", 
     col="firebrick", 
     xlab = "Time", 
     ylab = "Incidence per 100k")

plot(t , inc, type = "l", 
     col="firebrick", 
     xlab = "Time", 
     ylab = "Incidence per 100k",
     main = "Simulated TB incidence vs Paraguay estimate",
     xlim = c(950,1000),
     ylim = c(0,175))
points(999,62,pch=19)


plot(t , inc_p, type = "l", 
     col="navy", 
     xlab = "Time", 
     ylab = "Incidence per 100k",
     main = "Simulated TB incidence  in prisons vs Paraguay estimate")


plot(t , inc_p, type = "l", 
     col="navy", 
     xlab = "Time", 
     ylab = "Incidence per 100k",
     main = "Simulated TB incidence  in prisons vs Paraguay estimate",
     xlim = c(950,1000),
    ylim = c(0,8000))
points(999,4729,pch=19)

plot(t , prev_p, type = "l", 
     col="orange3", 
     xlab = "Time", 
     ylab = "Active TB prevalence (%)",
     main = "Simulated TB prevalence in prisons vs Cross-sectional estimate in Paraguay prisons",
     xlim = c(950,1000),
     ylim = c(0,15))
points(999,6.478,pch=19)

# Interventions -----------------------------------------------------------

# Baseline

pars<-list()

runs0<-update_intervention(sys,1001,init_state0,10,pars,"Baseline_100%")

runs0$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs0$incidence_community)) *1e5 

runs0$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs0$incidence_prison))/runs0$prison_pop*1e5 

runs0$cases_averted<-0

runs0$number_needed_screen<-0


# Algorithm 7 with 100% coverage
pars1=list(
  screen_exit_a = 0.818,
  screen_entry_a = 0.818,
  screen_yearly_a = 0.818,
  screen_exit_s = 0.909,
  screen_entry_s = 0.909,
  screen_yearly_s = 0.909
)

runs1<-update_intervention(sys,1001,init_state0,10,pars1,"Algorithm7_100%")

runs1$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs1$incidence_community)) *1e5 

runs1$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs1$incidence_prison))/runs0$prison_pop*1e5 

runs1$cases_averted<-runs0$n_cummulative_tb_cases- runs1$n_cummulative_tb_cases

runs1$number_needed_screen<-runs1$n_cummulative_screened/runs1$n_cumulative_detected

# Algorithm 7 with 75% coverage
pars2=list(
  screen_exit_a = 0.818*0.75,
  screen_entry_a = 0.818*0.75,
  screen_yearly_a = 0.818*0.75,
  screen_exit_s = 0.909*0.75,
  screen_entry_s = 0.909*0.75,
  screen_yearly_s = 0.909*0.75
)

runs2<-update_intervention(sys,1001,init_state0,10,pars2,"Algorithm7_75%")

runs2$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs2$incidence_community)) *1e5 

runs2$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs2$incidence_prison))/runs0$prison_pop*1e5 

runs2$cases_averted<-runs0$n_cummulative_tb_cases- runs2$n_cummulative_tb_cases

runs2$number_needed_screen<-runs2$n_cummulative_screened*0.75/runs2$n_cumulative_detected

# Algorithm 7 with 50% coverage
pars3=list(
  screen_exit_a = 0.818*0.5,
  screen_entry_a = 0.818*0.5,
  screen_yearly_a = 0.818*0.5,
  screen_exit_s = 0.909*0.5,
  screen_entry_s = 0.909*0.5,
  screen_yearly_s = 0.909*0.5
)

runs3<-update_intervention(sys,1001,init_state0,10,pars3,"Algorithm7_50%")

runs3$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs3$incidence_community)) *1e5 

runs3$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs3$incidence_prison))/runs0$prison_pop*1e5 

runs3$cases_averted<-runs0$n_cummulative_tb_cases- runs3$n_cummulative_tb_cases

runs3$number_needed_screen<-runs3$n_cummulative_screened*0.5/runs3$n_cumulative_detected

# Algorithm 15 with 100% coverage
pars4=list(
  screen_exit_a = 0.364,
  screen_entry_a = 0.364,
  screen_yearly_a = 0.364,
  screen_exit_s = 0.909,
  screen_entry_s = 0.909,
  screen_yearly_s = 0.909
)

runs4<-update_intervention(sys,1001,init_state0,10,pars4,"Algorithm15_100%")

runs4$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs4$incidence_community)) *1e5 

runs4$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs4$incidence_prison))/runs0$prison_pop*1e5 

runs4$cases_averted<-runs0$n_cummulative_tb_cases- runs4$n_cummulative_tb_cases

runs4$number_needed_screen<-runs4$n_cummulative_screened/runs4$n_cumulative_detected

# Algorithm 35 with 100% coverage
pars5=list(
  screen_exit_a = 0.818,
  screen_entry_a = 0.818,
  screen_yearly_a = 0.818,
  screen_exit_s = 0.800,
  screen_entry_s = 0.800,
  screen_yearly_s = 0.800
)

runs5<-update_intervention(sys,1001,init_state0,10,pars5,"Algorithm35_100%")

runs5$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs5$incidence_community)) *1e5 

runs5$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs5$incidence_prison))/runs0$prison_pop*1e5 

runs5$cases_averted<-runs0$n_cummulative_tb_cases- runs5$n_cummulative_tb_cases

runs5$number_needed_screen<-runs5$n_cummulative_screened/runs5$n_cumulative_detected

# Algorithm 37 with 100% coverage
pars6=list(
  screen_exit_a = 0.364,
  screen_entry_a = 0.364,
  screen_yearly_a = 0.364,
  screen_exit_s = 0.936,
  screen_entry_s = 0.936,
  screen_yearly_s = 0.936
)

runs6<-update_intervention(sys,1001,init_state0,10,pars6,"Algorithm37_100%")

runs6$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs6$incidence_community)) *1e5 

runs6$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs6$incidence_prison))/runs0$prison_pop*1e5 

runs6$cases_averted<-runs0$n_cummulative_tb_cases- runs6$n_cummulative_tb_cases

runs6$number_needed_screen<-runs6$n_cummulative_screened/runs6$n_cumulative_detected

# Algorithm 35 with 75% coverage
pars7=list(
  screen_exit_a = 0.818*0.75,
  screen_entry_a = 0.818*0.75,
  screen_yearly_a = 0.818*0.75,
  screen_exit_s = 0.800*0.75,
  screen_entry_s = 0.800*0.75,
  screen_yearly_s = 0.800*0.75
)

runs7<-update_intervention(sys,1001,init_state0,10,pars7,"Algorithm35_75%")

runs7$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs7$incidence_community)) *1e5 

runs7$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs7$incidence_prison))/runs0$prison_pop*1e5 

runs7$cases_averted<-runs0$n_cummulative_tb_cases- runs7$n_cummulative_tb_cases

runs7$number_needed_screen<-runs7$n_cummulative_screened*0.75/runs7$n_cumulative_detected

#algorithm 15 with 75% coverage
pars8=list(
  screen_exit_a = 0.364*0.75,
  screen_entry_a = 0.364*0.75,
  screen_yearly_a = 0.364*0.75,
  screen_exit_s = 0.909*0.75,
  screen_entry_s = 0.909*0.75,
  screen_yearly_s = 0.909*0.75
)

runs8<-update_intervention(sys,1001,init_state0,10,pars8,"Algorithm15_75%")

runs8$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs8$incidence_community)) *1e5 

runs8$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs8$incidence_prison))/runs0$prison_pop*1e5 

runs8$cases_averted<-runs0$n_cummulative_tb_cases- runs8$n_cummulative_tb_cases

runs8$number_needed_screen<-runs8$n_cummulative_screened*0.75/runs8$n_cumulative_detected

#algorithm 37 with 75% coverage
pars9=list(
  screen_exit_a = 0.364*0.75,
  screen_entry_a = 0.364*0.75,
  screen_yearly_a = 0.364*0.75,
  screen_exit_s = 0.936*0.75,
  screen_entry_s = 0.936*0.75,
  screen_yearly_s = 0.936*0.75
)

runs9<-update_intervention(sys,1001,init_state0,10,pars9,"Algorithm37_75%")

runs9$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs9$incidence_community)) *1e5 

runs9$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs9$incidence_prison))/runs0$prison_pop*1e5 

runs9$cases_averted<-runs0$n_cummulative_tb_cases- runs9$n_cummulative_tb_cases

runs9$number_needed_screen<-runs9$n_cummulative_screened*0.75/runs9$n_cumulative_detected

#algorithm 15 with 50% coverage
pars10=list(
  screen_exit_a = 0.364*0.5,
  screen_entry_a = 0.364*0.5,
  screen_yearly_a = 0.364*0.5,
  screen_exit_s = 0.909*0.5,
  screen_entry_s = 0.909*0.5,
  screen_yearly_s = 0.909*0.5
)

runs10<-update_intervention(sys,1001,init_state0,10,pars10,"Algorithm15_50%")

runs10$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs10$incidence_community)) *1e5 

runs10$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs10$incidence_prison))/runs0$prison_pop*1e5 

runs10$cases_averted<-runs0$n_cummulative_tb_cases- runs10$n_cummulative_tb_cases

runs10$number_needed_screen<-runs10$n_cummulative_screened*0.5/runs10$n_cumulative_detected

#algorithm 35 with 50% coverage
pars11=list(
  screen_exit_a = 0.818*0.5,
  screen_entry_a = 0.818*0.5,
  screen_yearly_a = 0.818*0.5,
  screen_exit_s = 0.800*0.5,
  screen_entry_s = 0.800*0.5,
  screen_yearly_s = 0.800*0.5
)

runs11<-update_intervention(sys,1001,init_state0,10,pars11,"Algorithm35_50%")

runs11$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs11$incidence_community)) *1e5 

runs11$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs11$incidence_prison))/runs0$prison_pop*1e5 

runs11$cases_averted<-runs0$n_cummulative_tb_cases- runs11$n_cummulative_tb_cases

runs11$number_needed_screen<-runs11$n_cummulative_screened*0.5/runs11$n_cumulative_detected

#algorithm 37 with 50% coverage
pars12=list(
  screen_exit_a = 0.364*0.5,
  screen_entry_a = 0.364*0.5,
  screen_yearly_a = 0.364*0.5,
  screen_exit_s = 0.936*0.5,
  screen_entry_s = 0.936*0.5,
  screen_yearly_s = 0.936*0.5
)

runs12<-update_intervention(sys,1001,init_state0,10,pars12,"Algorithm37_50%")

runs12$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs12$incidence_community)) *1e5 

runs12$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs12$incidence_prison))/runs0$prison_pop*1e5 

runs12$cases_averted<-runs0$n_cummulative_tb_cases- runs12$n_cummulative_tb_cases

runs12$number_needed_screen<-runs12$n_cummulative_screened*0.5/runs12$n_cumulative_detected

# Algorithm paraguay
pars13=list(
  screen_exit_a = 0.818 * 0.093,
  screen_entry_a = 0.818 * 0.093,
  screen_yearly_a = 0.818 * 0.093,
  screen_exit_s = 0.936 * 0.093,
  screen_entry_s = 0.936 * 0.093,
  screen_yearly_s = 0.936 * 0.093
)

runs13<-update_intervention(sys,1001,init_state0,10,pars13,"Algorithm paraguay_100%")

runs13$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs13$incidence_community)) *1e5 

runs13$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs13$incidence_prison))/runs0$prison_pop*1e5 

runs13$cases_averted<-runs0$n_cummulative_tb_cases- runs13$n_cummulative_tb_cases

runs13$number_needed_screen<-runs13$n_cummulative_screened*0.093/runs13$n_cumulative_detected

#algorithm 6 with 100% coverage
pars14=list(
  screen_exit_a = 0.818,
  screen_entry_a = 0.818,
  screen_yearly_a = 0.818,
  screen_exit_s = 0.936,
  screen_entry_s = 0.936,
  screen_yearly_s = 0.936
)

runs14<-update_intervention(sys,1001,init_state0,10,pars14,"Algorithm6_100%")

runs14$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs14$incidence_community)) *1e5 

runs14$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs14$incidence_prison))/runs0$prison_pop*1e5 

runs14$cases_averted<-runs0$n_cummulative_tb_cases- runs14$n_cummulative_tb_cases

runs14$number_needed_screen<-runs14$n_cummulative_screened/runs14$n_cumulative_detected

#algorithm 6 with 75% coverage
pars15=list(
  screen_exit_a = 0.818 * 0.75,
  screen_entry_a = 0.818 * 0.75,
  screen_yearly_a = 0.818 * 0.75,
  screen_exit_s = 0.936 * 0.75,
  screen_entry_s = 0.936 * 0.75,
  screen_yearly_s = 0.936 * 0.75
)

runs15<-update_intervention(sys,1001,init_state0,10,pars15,"Algorithm6_75%")

runs15$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs15$incidence_community)) *1e5 

runs15$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs15$incidence_prison))/runs0$prison_pop*1e5 

runs15$cases_averted<-runs0$n_cummulative_tb_cases- runs15$n_cummulative_tb_cases

runs15$number_needed_screen<-runs15$n_cummulative_screened/runs15$n_cumulative_detected

#algorithm 6 with 50% coverage
pars16=list(
  screen_exit_a = 0.818 * 0.5,
  screen_entry_a = 0.818 * 0.5,
  screen_yearly_a = 0.818 * 0.5,
  screen_exit_s = 0.936 * 0.5,
  screen_entry_s = 0.936 * 0.5,
  screen_yearly_s = 0.936 * 0.5
)

runs16<-update_intervention(sys,1001,init_state0,10,pars16,"Algorithm6_50%")

runs16$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs16$incidence_community)) *1e5 

runs16$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs16$incidence_prison))/runs0$prison_pop*1e5 

runs16$cases_averted<-runs0$n_cummulative_tb_cases- runs16$n_cummulative_tb_cases

runs16$number_needed_screen<-runs16$n_cummulative_screened/runs16$n_cumulative_detected

# Algorithm 24 with 100% coverage
pars17=list(
  screen_exit_a = 0.545,
  screen_entry_a = 0.545,
  screen_yearly_a = 0.545,
  screen_exit_s = 0.845,
  screen_entry_s = 0.845,
  screen_yearly_s = 0.845
)

runs17<-update_intervention(sys,1001,init_state0,10,pars17,"Algorithm24_100%")

runs17$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs17$incidence_community)) *1e5 

runs17$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs17$incidence_prison))/runs0$prison_pop*1e5 

runs17$cases_averted<-runs0$n_cummulative_tb_cases- runs17$n_cummulative_tb_cases

runs17$number_needed_screen<-runs17$n_cummulative_screened/runs17$n_cumulative_detected

# Algorithm 24 with 75% coverage
pars18=list(
  screen_exit_a = 0.545 * 0.75,
  screen_entry_a = 0.545 * 0.75,
  screen_yearly_a = 0.545 * 0.75,
  screen_exit_s = 0.845 * 0.75,
  screen_entry_s = 0.845 * 0.75,
  screen_yearly_s = 0.845 * 0.75
)

runs18<-update_intervention(sys,1001,init_state0,10,pars18,"Algorithm24_75%")

runs18$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs18$incidence_community)) *1e5 

runs18$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs18$incidence_prison))/runs0$prison_pop*1e5 

runs18$cases_averted<-runs0$n_cummulative_tb_cases- runs18$n_cummulative_tb_cases

runs18$number_needed_screen<-runs18$n_cummulative_screened/runs18$n_cumulative_detected

# Algorithm 24 with 50% coverage
pars19=list(
  screen_exit_a = 0.545 * 0.5,
  screen_entry_a = 0.545 * 0.5,
  screen_yearly_a = 0.545 * 0.5,
  screen_exit_s = 0.845 * 0.5,
  screen_entry_s = 0.845 * 0.5,
  screen_yearly_s = 0.845 * 0.5
)

runs19<-update_intervention(sys,1001,init_state0,10,pars19,"Algorithm24_50%")

runs19$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs19$incidence_community)) *1e5 

runs19$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs19$incidence_prison))/runs0$prison_pop*1e5 

runs19$cases_averted<-runs0$n_cummulative_tb_cases- runs19$n_cummulative_tb_cases

runs19$number_needed_screen<-runs19$n_cummulative_screened/runs19$n_cumulative_detected

# Algorithm 54 with 100% coverage
pars20=list(
  screen_exit_a = 0.545,
  screen_entry_a = 0.545,
  screen_yearly_a = 0.545,
  screen_exit_s = 0.864,
  screen_entry_s = 0.864,
  screen_yearly_s = 0.864
)

runs20<-update_intervention(sys,1001,init_state0,10,pars20,"Algorithm54_100%")

runs20$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs20$incidence_community)) *1e5 

runs20$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs20$incidence_prison))/runs0$prison_pop*1e5 

runs20$cases_averted<-runs0$n_cummulative_tb_cases- runs20$n_cummulative_tb_cases

runs20$number_needed_screen<-runs20$n_cummulative_screened/runs20$n_cumulative_detected

# Algorithm 54 with 75% coverage
pars21=list(
  screen_exit_a = 0.545 * 0.75,
  screen_entry_a = 0.545 * 0.75,
  screen_yearly_a = 0.545 * 0.75,
  screen_exit_s = 0.864 * 0.75,
  screen_entry_s = 0.864 * 0.75,
  screen_yearly_s = 0.864 * 0.75
)

runs21<-update_intervention(sys,1001,init_state0,10,pars21,"Algorithm54_75%")

runs21$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs21$incidence_community)) *1e5 

runs21$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs21$incidence_prison))/runs0$prison_pop*1e5 

runs21$cases_averted<-runs0$n_cummulative_tb_cases- runs21$n_cummulative_tb_cases

runs21$number_needed_screen<-runs21$n_cummulative_screened/runs21$n_cumulative_detected

# Algorithm 54 with 50% coverage
pars22=list(
  screen_exit_a = 0.545 * 0.5,
  screen_entry_a = 0.545 * 0.5,
  screen_yearly_a = 0.545 * 0.5,
  screen_exit_s = 0.864 * 0.5,
  screen_entry_s = 0.864 * 0.5,
  screen_yearly_s = 0.864 * 0.5
)

runs22<-update_intervention(sys,1001,init_state0,10,pars22,"Algorithm54_50%")

runs22$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs22$incidence_community)) *1e5 

runs22$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs22$incidence_prison))/runs0$prison_pop*1e5 

runs22$cases_averted<-runs0$n_cummulative_tb_cases- runs22$n_cummulative_tb_cases

runs22$number_needed_screen<-runs22$n_cummulative_screened/runs22$n_cumulative_detected

# Algorithm 56 with 100% coverage, annual screening = every 3 months (first BMI selection)
pars23=list(
  screen_exit_a = 0.364,
  screen_entry_a = 0.364,
  screen_yearly_a = 0.364 * 4,
  screen_exit_s = 0.627,
  screen_entry_s = 0.627,
  screen_yearly_s = 0.627 * 4
)

runs23<-update_intervention(sys,1001,init_state0,10,pars23,"Algorithm56_100%")

runs23$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs23$incidence_community)) *1e5 

runs23$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs23$incidence_prison))/runs0$prison_pop*1e5 

runs23$cases_averted<-runs0$n_cummulative_tb_cases- runs23$n_cummulative_tb_cases

runs23$number_needed_screen<-runs23$n_cummulative_screened/runs23$n_cumulative_detected

# Algorithm 56 with 75% coverage, annual screening = every 3 months (first BMI selection)
pars24=list(
  screen_exit_a = 0.364 * 0.75,
  screen_entry_a = 0.364 * 0.75,
  screen_yearly_a = 0.364 * 4 * 0.75,
  screen_exit_s = 0.627 * 0.75,
  screen_entry_s = 0.627 * 0.75,
  screen_yearly_s = 0.627 * 4 * 0.75
)

runs24<-update_intervention(sys,1001,init_state0,10,pars24,"Algorithm56_75%")

runs24$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs24$incidence_community)) *1e5 

runs24$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs24$incidence_prison))/runs0$prison_pop*1e5 

runs24$cases_averted<-runs0$n_cummulative_tb_cases- runs24$n_cummulative_tb_cases

runs24$number_needed_screen<-runs24$n_cummulative_screened/runs24$n_cumulative_detected

# Algorithm 56 with 75% coverage, annual screening = every 3 months (first BMI selection)
pars25=list(
  screen_exit_a = 0.364 * 0.5,
  screen_entry_a = 0.364 * 0.5,
  screen_yearly_a = 0.364 * 4 * 0.5,
  screen_exit_s = 0.627 * 0.5,
  screen_entry_s = 0.627 * 0.5,
  screen_yearly_s = 0.627 * 4 * 0.5
)

runs25<-update_intervention(sys,1001,init_state0,10,pars25,"Algorithm56_50%")

runs25$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs25$incidence_community)) *1e5 

runs25$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs25$incidence_prison))/runs0$prison_pop*1e5 

runs25$cases_averted<-runs0$n_cummulative_tb_cases- runs25$n_cummulative_tb_cases

runs25$number_needed_screen<-runs25$n_cummulative_screened/runs25$n_cumulative_detected

# Algorithm 57 with 100% coverage, annual screening = every 3 months (first BMI selection)
pars26=list(
  screen_exit_a = 0.364,
  screen_entry_a = 0.364,
  screen_yearly_a = 0.364 * 4,
  screen_exit_s = 0.636,
  screen_entry_s = 0.636,
  screen_yearly_s = 0.636 * 4
)

runs26<-update_intervention(sys,1001,init_state0,10,pars26,"Algorithm57_100%")

runs26$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs26$incidence_community)) *1e5 

runs26$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs26$incidence_prison))/runs0$prison_pop*1e5 

runs26$cases_averted<-runs0$n_cummulative_tb_cases- runs26$n_cummulative_tb_cases

runs26$number_needed_screen<-runs26$n_cummulative_screened/runs26$n_cumulative_detected

# Algorithm 57 with 75% coverage, annual screening = every 3 months (first BMI selection)
pars27=list(
  screen_exit_a = 0.364 * 0.75,
  screen_entry_a = 0.364* 0.75,
  screen_yearly_a = 0.364 * 4 * 0.75,
  screen_exit_s = 0.636 * 0.75,
  screen_entry_s = 0.636 * 0.75,
  screen_yearly_s = 0.636 * 4 * 0.75
)

runs27<-update_intervention(sys,1001,init_state0,10,pars27,"Algorithm57_75%")

runs27$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs27$incidence_community)) *1e5 

runs27$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs27$incidence_prison))/runs0$prison_pop*1e5 

runs27$cases_averted<-runs0$n_cummulative_tb_cases- runs27$n_cummulative_tb_cases

runs27$number_needed_screen<-runs27$n_cummulative_screened/runs27$n_cumulative_detected

# Algorithm 57 with 50% coverage, annual screening = every 3 months (first BMI selection)
pars28=list(
  screen_exit_a = 0.364 * 0.5,
  screen_entry_a = 0.364* 0.5,
  screen_yearly_a = 0.364 * 4 * 0.5,
  screen_exit_s = 0.636 * 0.5,
  screen_entry_s = 0.636 * 0.5,
  screen_yearly_s = 0.636 * 4 * 0.5
)

runs28<-update_intervention(sys,1001,init_state0,10,pars28,"Algorithm57_50%")

runs28$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs28$incidence_community)) *1e5 

runs28$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs28$incidence_prison))/runs0$prison_pop*1e5 

runs28$cases_averted<-runs0$n_cummulative_tb_cases- runs28$n_cummulative_tb_cases

runs28$number_needed_screen<-runs28$n_cummulative_screened/runs28$n_cumulative_detected

# Join together all results 
df<-rbind(runs0,runs1,runs2,runs3, runs4, runs5, runs6, runs7,runs8,runs9,runs10,runs11,runs12,runs14, runs15, runs16,
          runs17, runs18, runs19, runs20, runs21, runs22, runs26, runs24, runs25, runs26, runs27, runs28)
df<- df %>%
  separate(scenario, into = c("Scenario", "Coverage"), sep = "_")

#Plot TB community 
ggplot(df, aes(x=years,y=incidence_community_rate,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "TB incidence in the community", y="Incidence per 100,000")+
  ylim(30,65)+
  theme_minimal()

#Plot TB prison 
ggplot(df, aes(x=years,y=incidence_prison_rate,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "TB incidence in the prison", y="Incidence per 100,000")+
  ylim(0,5.5e3)+
  theme_minimal()


ggplot(df, aes(x=years,y=n_annual_detected,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "Annual TB cases detected", y="TB cases")+
  ylim(0,750)+
  theme_minimal()


ggplot(df, aes(x=years,y=cases_averted,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "Cumulative TB cases averted", y="TB cases")+
  ylim(0,15000)+
  theme_minimal()


ggplot(df, aes(x=years,y=number_needed_screen,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "Number needed to screen", y="Screened/detected")+
  ylim(0, 200)+
  xlim(1,10)+
  theme_minimal()




