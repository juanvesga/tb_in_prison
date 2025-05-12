set.seed(42)

## Load packages

library(odin2)
library(dust2)
library(monty)


# deterministic model general population tb
uli_ode <- odin({
  

# Community model ---------------------------------------------------------

  deriv(U)  <- 
    births - 
    U * lambda - 
    U * mu - 
    U * r_incar + 
    U_p * r_release
  
  deriv(L)  <- 
    U * lambda * (1-fast) + 
    R * (lambda * (1-fast) * imm) - 
    L * (mu + slow) - 
    L * r_incar + 
    L_p * r_release
  
  deriv(Ia) <- 
    U * lambda * fast + 
    R * (lambda * fast * imm) + 
    L * slow - 
    Ia * (sigma + mu) - 
    Ia * r_incar + 
    Ia_p * r_release
  
  deriv(Is) <- 
    Ia * sigma - 
    Is * (mutb + mu + self_cure + r_tx) - 
    Is * r_incar + 
    Is_p * r_release
  
  deriv(R)  <- 
    Is*(self_cure + r_tx) - 
    R * (imm*lambda + mu) - 
    R * r_incar + 
    R_p * r_release
  
  deriv(incidence) <-  
    U * lambda * fast + 
    R * lambda * fast * imm +  
    L * slow 
  

# Prison model ------------------------------------------------------------
  
  deriv(U_p)  <- 
    - U_p * (lambda_p + U_p) +
    U   * r_incar - 
    U_p * r_release
  
  deriv(L_p)  <- 
    U_p * lambda_p * (1-fast_p) + 
    R_p * (lambda_p * (1-fast_p) * imm) - 
    L_p * (mu + slow_p) +
    L * r_incar -
    L_p * r_release
  
  deriv(Ia_p) <- 
    U_p * lambda_p * fast_p + 
    R_p * (lambda_p * fast_p * imm) +  
    L_p * slow_p - Ia_p * (sigma_p + mu) + 
    Ia * r_incar - 
    Ia_p * r_release
  
  deriv(Is_p) <- 
    Ia_p * sigma_p - 
    Is_p * (mutb + mu + self_cure + r_tx_p) + 
    Is * r_incar - 
    Is_p * r_release
  
  deriv(R_p)  <- 
    Is_p * (self_cure + r_tx_p) - 
    R_p * (imm*lambda_p + mu) + 
    R * r_incar - 
    R_p * r_release
  
  deriv(incidence_p) <-  
    U_p * lambda_p * fast_p + 
    R_p * lambda_p * fast_p * imm +  
    L_p * slow_p 
  
  N       <- U + L + Ia + Is + R
  N_p     <- U_p + L_p + Ia_p + Is_p + R_p
  births  <- mu*(N+N_p) + mutb*(Is+Is_p) 
  lambda  <- beta   * (Ia  +Is  )/N # force of infection
  lambda_p<- beta_p * (Ia_p+Is_p)/N_p # force of infection prisons
  
  # Known Model Parameters
  l_exp    <- 72              # Life expectancy
  tb_dur    <- 3               # Duration of infectious period (years)
  mu       <- 1/l_exp         # Background mortality rate
  mutb     <- 0.5*(1/tb_dur)  # TB mortality rate
  self_cure<- 0.5*(1/tb_dur)  # recovery
  fast     <- 0.1             # Fraction fast progressing to active TB
  slow     <- 0.05*1/l_exp    # Remote reactivation
  sigma    <- 1/0.5           # symptom development (6 mo)
  fast_p   <- 0.18             # Fraction fast progressing to active TB
  slow_p   <- 0.08*1/l_exp    # Remote reactivation
  sigma_p  <- 1/0.33          # symptom development (4 mo)
  imm      <- 0.5             # Infectiousness decline (partial immunity)
  r_tx     <- if (time > 950) 1*0.7*0.7 else 0#       # Careseeking rate (1 year)
  r_tx_p   <- if (time > 950) prison_tx  else 0#       # Careseeking rate (1 year)
  I0       <- 1e-6
  P0       <- 250/1e5
  r_incar  <- 0.001          # rate of incarceration
  r_release<- 1/2.5          # rate of release (1/mean prison term)


# External inputs ---------------------------------------------------------
  beta     <-parameter(5)
  beta_p   <-parameter(5)
  prison_tx <- parameter(0.49)

  
    
  initial(U) <- 1-I0-P0
  initial(L) <- 0
  initial(Ia) <- 0
  initial(Is) <- I0
  initial(R) <- 0
  
  initial(U_p) <- P0-(I0*P0)
  initial(L_p) <- 0
  initial(Ia_p) <- 0
  initial(Is_p) <- I0*P0
  initial(R_p) <- 0
  
  initial(incidence) <- 0 
  initial(incidence_p) <- 0 
  
})



sys <- dust_system_create(uli_ode, pars = list(beta=4.5,beta_p=34.5))

dust_system_set_state_initial(sys)
t <- seq(0, 1000)
y <- dust_system_simulate(sys, t)

y <- dust_unpack_state(sys, y)

totalpop <- y$Is + y$Ia + y$U + y$L + y$R + y$Is_p + y$Ia_p + y$U_p + y$L_p + y$R_p
prison_pop<-y$Is_p + y$Ia_p + y$U_p + y$L_p + y$R_p

inc<-c(0,diff(y$incidence)*1e5)

inc_p<-c(0,diff(y$incidence_p))/prison_pop*1e5



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
points(999,46,pch=19)


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




