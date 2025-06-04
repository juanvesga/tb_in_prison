set.seed(42)

## Load packages

library(odin2)
library(dust2)
library(monty)
library(here)
library(ggplot2)

root<-here()

# Call model script
source(file.path(root,"v2 model.R"))


# Create model object
sys <- dust_system_create(uli_ode, 
                          pars = list(beta=3,beta_p=16))

# Set default initial state of the model 
dust_system_set_state_initial(sys)

# Run simulation from time zero
t <- seq(0, 1000)
y <- dust_system_simulate(sys, t)

# Get model state at the latest point for future simulations
init_state0<-dust_system_state(sys)

y <- dust_unpack_state(sys, y)

names(init_state0)<-paste0(names(y))

totalpop <- y$Is + y$Ia + y$U + y$L + y$R + y$Is_p + y$Ia_p + y$U_p + y$L_p + y$R_p + y$Tr + y$Tr_p

prison_pop<-y$Is_p + y$Ia_p + y$U_p + y$L_p + y$R_p + y$Tr_p

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
points(999,4200,pch=19)



# Interventions -----------------------------------------------------------

# Baseline

pars<-list()

runs0<-update_intervention(sys,1001,init_state0,10,pars,"Baseline")

runs0$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs0$incidence_community)) *1e5 

runs0$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs0$incidence_prison))/runs0$prison_pop*1e5 

runs0$cases_averted<-0

runs0$number_needed_screen<-0


# Algorithm 1
pars1=list(
  screen_exit_a = 0.818,
  screen_entry_a = 0.818,
  screen_yearly_a = 0.818,
  screen_exit_s = 0.891,
  screen_entry_s = 0.891,
  screen_yearly_s = 0.891
)

runs1<-update_intervention(sys,1001,init_state0,10,pars1,"Algorithm 1")

runs1$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs1$incidence_community)) *1e5 

runs1$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs1$incidence_prison))/runs0$prison_pop*1e5 

runs1$cases_averted<-runs0$n_cummulative_tb_cases- runs1$n_cummulative_tb_cases

runs1$number_needed_screen<-runs1$n_cummulative_screened/runs1$n_cumulative_detected

# Algorithm 2
pars2=list(
  screen_exit_a = 0.364,
  screen_entry_a = 0.364,
  screen_yearly_a = 0.364,
  screen_exit_s = 0.882,
  screen_entry_s = 0.882,
  screen_yearly_s = 0.882
)

runs2<-update_intervention(sys,1001,init_state0,10,pars2,"Algorithm 2")

runs2$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs2$incidence_community)) *1e5 

runs2$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs2$incidence_prison))/runs0$prison_pop*1e5 

runs2$cases_averted<-runs0$n_cummulative_tb_cases- runs2$n_cummulative_tb_cases

runs2$number_needed_screen<-runs2$n_cummulative_screened/runs2$n_cumulative_detected

# Algorithm 3
pars3=list(
  screen_exit_a = 0.818,
  screen_entry_a = 0.818,
  screen_yearly_a = 0.818,
  screen_exit_s = 0.765,
  screen_entry_s = 0.765,
  screen_yearly_s = 0.765
)

runs3<-update_intervention(sys,1001,init_state0,10,pars3,"Algorithm 3")

runs3$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs3$incidence_community)) *1e5 

runs3$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs3$incidence_prison))/runs0$prison_pop*1e5 

runs3$cases_averted<-runs0$n_cummulative_tb_cases- runs3$n_cummulative_tb_cases

runs3$number_needed_screen<-runs3$n_cummulative_screened/runs3$n_cumulative_detected

# Algorithm 4
pars4=list(
  screen_exit_a = 0.364,
  screen_entry_a = 0.364,
  screen_yearly_a = 0.364,
  screen_exit_s = 0.924,
  screen_entry_s = 0.924,
  screen_yearly_s = 0.924
)

runs4<-update_intervention(sys,1001,init_state0,10,pars4,"Algorithm 4")

runs4$incidence_community_rate<-diff( c(y$incidence[length(y$incidence)-1],runs4$incidence_community)) *1e5 

runs4$incidence_prison_rate<-diff( c(y$incidence_p[length(y$incidence_p)-1],runs4$incidence_prison))/runs0$prison_pop*1e5 

runs4$cases_averted<-runs0$n_cummulative_tb_cases- runs4$n_cummulative_tb_cases

runs4$number_needed_screen<-runs4$n_cummulative_screened/runs4$n_cumulative_detected


# Join together all results 
df<-rbind(runs0,runs1,runs2, runs3, runs4)

#Plot TB community 
ggplot(df, aes(x=years,y=incidence_community_rate,colour = scenario))+
  geom_line()+
  labs(title = "TB incidence in the community", y="Incidence per 100,000")+
  ylim(0,75)+
  theme_minimal()

#Plot TB prison 
ggplot(df, aes(x=years,y=incidence_prison_rate,colour = scenario))+
  geom_line()+
  labs(title = "TB incidence in the prison", y="Incidence per 100,000")+
  ylim(0,5e3)+
  theme_minimal()


ggplot(df, aes(x=years,y=n_annual_detected,colour = scenario))+
  geom_line()+
  labs(title = "Annual TB cases detected", y="TB cases")+
  ylim(0,800)+
  theme_minimal()


ggplot(df, aes(x=years,y=cases_averted,colour = scenario))+
  geom_line()+
  labs(title = "Cumulative TB cases averted", y="TB cases")+
  ylim(0,20000)+
  theme_minimal()


ggplot(df, aes(x=years,y=number_needed_screen,colour = scenario))+
  geom_line()+
  labs(title = "Number needed to screen", y="Screened/detected")+
  ylim(0,200)+
  xlim(1,10)+
  theme_minimal()



