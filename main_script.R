set.seed(42)

## Load packages

library(odin2)
library(dust2)
library(monty)


# deterministic model general population tb
uli_ode <- odin({
  deriv(U)  <- births - U * lambda - mu*U 
  deriv(L)  <- U * lambda * (1-fast) + R * (lambda * (1-fast) * imm) - L * (mu + slow)
  deriv(Ia) <- U * lambda * fast + R * (lambda * fast * imm) +  L * slow - Ia * (sigma + mu)
  deriv(Is) <- Ia * sigma - Is * (mutb + mu + self_cure + r_tx)
  deriv(R)  <- Is*(self_cure + r_tx) - R * (imm*lambda + mu)
  deriv(incidence) <-  U * lambda * fast + R * lambda * fast * imm +  L * slow
  
  N     <- U+L+Ia+Is+R
  births<- mu*N + mutb*Is 
  lambda<- beta *(Ia+Is)/N # force of infection
  
  # Known Model Parameters
  l_exp    <- 72              # Life expectancy
  tb_dur    <- 3               # Duration of infectious period (years)
  mu       <- 1/l_exp         # Background mortality rate
  mutb     <- 0.5*(1/tb_dur)  # TB mortality rate
  self_cure<- 0.5*(1/tb_dur)  # recovery
  fast     <- 0.1             # Fraction fast progressing to active TB
  slow     <- 0.05*1/l_exp    # Remote reactivation
  sigma    <- 1/0.5           # symptom development (6 mo)
  imm      <- 0.5             # Infectiousness decline (partial immunity)
  r_tx     <- if (time > 950) 1*0.8*0.9 else 0#       # Careseeking rate (1 year)
  I0       <- 1e-6
  beta     <-parameter(5)
  
  initial(U) <- 1-I0
  initial(L) <- 0
  initial(Ia) <- 0
  initial(Is) <- I0
  initial(R) <- 0
  initial(incidence) <- 0 
  
})


sys <- dust_system_create(uli_ode, pars = list())

dust_system_set_state_initial(sys)
t <- seq(0, 1000)
y <- dust_system_simulate(sys, t)

dim(y)
y <- dust_unpack_state(sys, y)

inc<-c(0,diff(y$incidence)*1e5)

totalpop <- y$Is + y$Ia + y$U + y$L + y$R

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
