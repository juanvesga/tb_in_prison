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
     ylim = c(0,200))
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
     xlim = c(950,1000))
#ylim = c(0,200))
points(999,4729,pch=19)



# Interventions -----------------------------------------------------------

# Baseline

pars<-list()

runs0<-update_intervention(sys,1001,init_state0,10,pars,"Baseline_100%")

runs0$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs0$incidence_community)) *1e5 

runs0$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs0$incidence_prison))/runs0$prison_pop*1e5 

runs0$cases_averted<-0

runs0$number_needed_screen<-0


# Algorithm 1 with 100% coverage
pars1=list(
  screen_exit_a = 0.818,
  screen_entry_a = 0.818,
  screen_yearly_a = 0.818,
  screen_exit_s = 0.891,
  screen_entry_s = 0.891,
  screen_yearly_s = 0.891
)

runs1<-update_intervention(sys,1001,init_state0,10,pars1,"Algorithm1_100%")

runs1$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs1$incidence_community)) *1e5 

runs1$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs1$incidence_prison))/runs0$prison_pop*1e5 

runs1$cases_averted<-runs0$n_cummulative_tb_cases- runs1$n_cummulative_tb_cases

runs1$number_needed_screen<-runs1$n_cummulative_screened/runs1$n_cumulative_detected

# Algorithm 1 with 75% coverage
pars2=list(
  screen_exit_a = 0.818*0.75,
  screen_entry_a = 0.818*0.75,
  screen_yearly_a = 0.818*0.75,
  screen_exit_s = 0.891*0.75,
  screen_entry_s = 0.891*0.75,
  screen_yearly_s = 0.891*0.75
)

runs2<-update_intervention(sys,1001,init_state0,10,pars2,"Algorithm1_75%")

runs2$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs2$incidence_community)) *1e5 

runs2$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs2$incidence_prison))/runs0$prison_pop*1e5 

runs2$cases_averted<-runs0$n_cummulative_tb_cases- runs2$n_cummulative_tb_cases

runs2$number_needed_screen<-runs2$n_cummulative_screened/runs2$n_cumulative_detected

# Algorithm 1 with 50% coverage
pars3=list(
  screen_exit_a = 0.818*0.5,
  screen_entry_a = 0.818*0.5,
  screen_yearly_a = 0.818*0.5,
  screen_exit_s = 0.891*0.5,
  screen_entry_s = 0.891*0.5,
  screen_yearly_s = 0.891*0.5
)

runs3<-update_intervention(sys,1001,init_state0,10,pars3,"Algorithm1_50%")

runs3$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs3$incidence_community)) *1e5 

runs3$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs3$incidence_prison))/runs0$prison_pop*1e5 

runs3$cases_averted<-runs0$n_cummulative_tb_cases- runs3$n_cummulative_tb_cases

runs3$number_needed_screen<-runs3$n_cummulative_screened/runs3$n_cumulative_detected

# Algorithm 2 with 100% coverage
pars4=list(
  screen_exit_a = 0.364,
  screen_entry_a = 0.364,
  screen_yearly_a = 0.364,
  screen_exit_s = 0.882,
  screen_entry_s = 0.882,
  screen_yearly_s = 0.882
)

runs4<-update_intervention(sys,1001,init_state0,10,pars4,"Algorithm2_100%")

runs4$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs4$incidence_community)) *1e5 

runs4$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs4$incidence_prison))/runs0$prison_pop*1e5 

runs4$cases_averted<-runs0$n_cummulative_tb_cases- runs4$n_cummulative_tb_cases

runs4$number_needed_screen<-runs4$n_cummulative_screened/runs4$n_cumulative_detected

# Algorithm 3 with 100% coverage
pars5=list(
  screen_exit_a = 0.818,
  screen_entry_a = 0.818,
  screen_yearly_a = 0.818,
  screen_exit_s = 0.765,
  screen_entry_s = 0.765,
  screen_yearly_s = 0.765
)

runs5<-update_intervention(sys,1001,init_state0,10,pars5,"Algorithm3_100%")

runs5$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs5$incidence_community)) *1e5 

runs5$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs5$incidence_prison))/runs0$prison_pop*1e5 

runs5$cases_averted<-runs0$n_cummulative_tb_cases- runs5$n_cummulative_tb_cases

runs5$number_needed_screen<-runs5$n_cummulative_screened/runs5$n_cumulative_detected

# Algorithm 4 with 100% coverage
pars6=list(
  screen_exit_a = 0.364,
  screen_entry_a = 0.364,
  screen_yearly_a = 0.364,
  screen_exit_s = 0.924,
  screen_entry_s = 0.924,
  screen_yearly_s = 0.924
)

runs6<-update_intervention(sys,1001,init_state0,10,pars6,"Algorithm4_100%")

runs6$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs6$incidence_community)) *1e5 

runs6$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs6$incidence_prison))/runs0$prison_pop*1e5 

runs6$cases_averted<-runs0$n_cummulative_tb_cases- runs6$n_cummulative_tb_cases

runs6$number_needed_screen<-runs6$n_cummulative_screened/runs6$n_cumulative_detected

# Algorithm 3 with 75% coverage
pars7=list(
  screen_exit_a = 0.818*0.75,
  screen_entry_a = 0.818*0.75,
  screen_yearly_a = 0.818*0.75,
  screen_exit_s = 0.765*0.75,
  screen_entry_s = 0.765*0.75,
  screen_yearly_s = 0.765*0.75
)

runs7<-update_intervention(sys,1001,init_state0,10,pars7,"Algorithm3_75%")

runs7$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs7$incidence_community)) *1e5 

runs7$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs7$incidence_prison))/runs0$prison_pop*1e5 

runs7$cases_averted<-runs0$n_cummulative_tb_cases- runs7$n_cummulative_tb_cases

runs7$number_needed_screen<-runs7$n_cummulative_screened/runs7$n_cumulative_detected

#algorithm 2 with 75% coverage
pars8=list(
  screen_exit_a = 0.364*0.75,
  screen_entry_a = 0.364*0.75,
  screen_yearly_a = 0.364*0.75,
  screen_exit_s = 0.882*0.75,
  screen_entry_s = 0.882*0.75,
  screen_yearly_s = 0.882*0.75
)

runs8<-update_intervention(sys,1001,init_state0,10,pars8,"Algorithm2_75%")

runs8$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs8$incidence_community)) *1e5 

runs8$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs8$incidence_prison))/runs0$prison_pop*1e5 

runs8$cases_averted<-runs0$n_cummulative_tb_cases- runs8$n_cummulative_tb_cases

runs8$number_needed_screen<-runs8$n_cummulative_screened/runs8$n_cumulative_detected

#algorithm 4 with 75% coverage
pars9=list(
  screen_exit_a = 0.364*0.75,
  screen_entry_a = 0.364*0.75,
  screen_yearly_a = 0.364*0.75,
  screen_exit_s = 0.924*0.75,
  screen_entry_s = 0.924*0.75,
  screen_yearly_s = 0.924*0.75
)

runs9<-update_intervention(sys,1001,init_state0,10,pars9,"Algorithm4_75%")

runs9$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs9$incidence_community)) *1e5 

runs9$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs9$incidence_prison))/runs0$prison_pop*1e5 

runs9$cases_averted<-runs0$n_cummulative_tb_cases- runs9$n_cummulative_tb_cases

runs9$number_needed_screen<-runs9$n_cummulative_screened/runs9$n_cumulative_detected

#algorithm 2 with 50% coverage
pars10=list(
  screen_exit_a = 0.364*0.5,
  screen_entry_a = 0.364*0.5,
  screen_yearly_a = 0.364*0.5,
  screen_exit_s = 0.882*0.5,
  screen_entry_s = 0.882*0.5,
  screen_yearly_s = 0.882*0.5
)

runs10<-update_intervention(sys,1001,init_state0,10,pars10,"Algorithm2_50%")

runs10$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs10$incidence_community)) *1e5 

runs10$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs10$incidence_prison))/runs0$prison_pop*1e5 

runs10$cases_averted<-runs0$n_cummulative_tb_cases- runs10$n_cummulative_tb_cases

runs10$number_needed_screen<-runs10$n_cummulative_screened/runs10$n_cumulative_detected

#algorithm 3 with 50% coverage
pars11=list(
  screen_exit_a = 0.818*0.5,
  screen_entry_a = 0.818*0.5,
  screen_yearly_a = 0.818*0.5,
  screen_exit_s = 0.765*0.5,
  screen_entry_s = 0.765*0.5,
  screen_yearly_s = 0.765*0.5
)

runs11<-update_intervention(sys,1001,init_state0,10,pars11,"Algorithm3_50%")

runs11$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs11$incidence_community)) *1e5 

runs11$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs11$incidence_prison))/runs0$prison_pop*1e5 

runs11$cases_averted<-runs0$n_cummulative_tb_cases- runs11$n_cummulative_tb_cases

runs11$number_needed_screen<-runs11$n_cummulative_screened/runs11$n_cumulative_detected

#algorithm 4 with 50% coverage
pars12=list(
  screen_exit_a = 0.364*0.5,
  screen_entry_a = 0.364*0.5,
  screen_yearly_a = 0.364*0.5,
  screen_exit_s = 0.924*0.5,
  screen_entry_s = 0.924*0.5,
  screen_yearly_s = 0.924*0.5
)

runs12<-update_intervention(sys,1001,init_state0,10,pars12,"Algorithm4_50%")

runs12$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs12$incidence_community)) *1e5 

runs12$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs12$incidence_prison))/runs0$prison_pop*1e5 

runs12$cases_averted<-runs0$n_cummulative_tb_cases- runs12$n_cummulative_tb_cases

runs12$number_needed_screen<-runs12$n_cummulative_screened/runs12$n_cumulative_detected

# Join together all results 
df<-rbind(runs0,runs1,runs2,runs3, runs4, runs5, runs6, runs7,runs8,runs9,runs10,runs11,runs12)
df<- df %>%
  separate(scenario, into = c("Scenario", "Coverage"), sep = "_")

#Plot TB community 
ggplot(df, aes(x=years,y=incidence_community_rate,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "TB incidence in the community", y="Incidence per 100,000")+
  ylim(0,75)+
  theme_minimal()

#Plot TB prison 
ggplot(df, aes(x=years,y=incidence_prison_rate,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "TB incidence in the prison", y="Incidence per 100,000")+
  ylim(0,5e3)+
  theme_minimal()


ggplot(df, aes(x=years,y=n_annual_detected,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "Annual TB cases detected", y="TB cases")+
  ylim(0,800)+
  theme_minimal()


ggplot(df, aes(x=years,y=cases_averted,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "Cumulative TB cases averted", y="TB cases")+
  ylim(0,20000)+
  theme_minimal()


ggplot(df, aes(x=years,y=number_needed_screen,colour = Scenario, linetype = Coverage))+
  geom_line()+
  labs(title = "Number needed to screen", y="Screened/detected")+
  ylim(0,200)+
  xlim(1,10)+
  theme_minimal()




