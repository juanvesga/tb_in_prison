library(odin2)
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
    Ia_p * r_release *(1-screen_exit_a) +
    Tr * tx_duration * (1-tx_completion) * (1-tx_forgiveness)
  
  deriv(Is) <- 
    Ia * sigma - 
    Is * (mutb + mu + self_cure + (tx_initiation_delay*case_detection)) - 
    Is * r_incar + 
    Is_p * r_release *(1-screen_exit_s) 
  
  deriv(Tr) <-
    Is * tx_initiation_delay * case_detection -
    Tr * tx_duration +
    Is_p * r_release * screen_exit_s * screen_delay +
    Ia_p * r_release * screen_exit_a * screen_delay -
    Tr * mu -
    Tr * r_incar + 
    Tr_p * r_release 
  
  deriv(R)  <- 
    Is * (self_cure) - 
    R * (imm*lambda + mu) - 
    R * r_incar + 
    R_p * r_release +
    Tr * tx_duration * tx_completion +
    Tr * tx_duration * (1-tx_completion) * tx_forgiveness
  
  deriv(incidence) <-  
    U * lambda * fast + 
    R * lambda * fast * imm +  
    L * slow 
  
  
  # Prison model ------------------------------------------------------------
  
  deriv(U_p)  <- 
    - U_p * (lambda_p + mu) +
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
    L_p * slow_p - 
    Ia_p * (sigma_p + mu) + 
    Ia * r_incar * (1-screen_entry_a) - 
    Ia_p * r_release -
    Ia_p * screen_yearly_a * screen_delay +
    Tr_p * tx_duration * (1-tx_completion) * (1-tx_forgiveness)
  
  deriv(Is_p) <- 
    Ia_p * sigma_p - 
    Is_p * (mutb + mu + self_cure + (tx_iniation_delay_p * case_detection_p)) + 
    Is * r_incar * (1-screen_entry_s) - 
    Is_p * r_release -
    Is_p * screen_yearly_s * screen_delay
  
  deriv(Tr_p) <-
    Is_p * tx_iniation_delay_p * case_detection_p -
    Tr_p * mu -
    Tr_p * tx_duration +
    Tr * r_incar -
    Tr_p * r_release +
    Ia * r_incar * screen_entry_a * screen_delay +
    Is * r_incar * screen_entry_s * screen_delay +
    Ia_p * screen_yearly_a * screen_delay +
    Is_p * screen_yearly_s * screen_delay
  
  deriv(R_p)  <- 
    Is_p * self_cure - 
    R_p * (imm*lambda_p + mu) + 
    R * r_incar -
    R_p * r_release +
    Tr_p * tx_duration * tx_completion +
    Tr_p * tx_duration * (1-tx_completion) * tx_forgiveness
  
  
  # Output variables --------------------------------------------------------
  
  deriv(incidence_p) <-  
    U_p * lambda_p * fast_p + 
    R_p * lambda_p * fast_p * imm +  
    L_p * slow_p 
  
  deriv(new_imprisoned) <- N  * r_incar
  
  deriv(new_released)   <- N_p  * r_release
  
  deriv(new_detected)  <- 
    Ia * r_incar * screen_entry_a +
    Is * r_incar * screen_entry_s +
    Ia_p * screen_yearly_a +
    Is_p * screen_yearly_s +
    Is_p * r_release * screen_exit_s +
    Ia_p * r_release * screen_exit_a 
  
  
  # Population and force of infection ---------------------------------------
  
  N       <- U + L + Ia + Is + Tr + R
  
  N_p     <- U_p + L_p + Ia_p + Is_p + Tr_p + R_p
  
  births  <- mu*(N+N_p) + mutb*(Is+Is_p) 
  
  lambda  <- (1-mix_p*prob_inf_mix)*(beta   * (Ia  +Is  )/N) + 
    mix_p*prob_inf_mix*(beta_p * (Ia_p+Is_p)/N_p) # force of infection
  
  lambda_p<- (1-mix_p*prob_inf_mix)*(beta_p * (Ia_p+Is_p)/N_p) + 
    mix_p*prob_inf_mix*(beta * (Ia+Is)/N) # force of infection prisons
  
  
  # Model parameters (fixed) ------------------------------------------------
  
  l_exp    <- 70             # Life expectancy
  tb_dur    <- 3             # Duration of infectious period (years)
  mu       <- 1/l_exp        # Background mortality rate
  mutb     <- 0.5*(1/tb_dur) # TB mortality rate
  self_cure<- 0.5*(1/tb_dur) # recovery
  fast     <- 0.1            # Fraction fast progressing to active TB
  slow     <- 0.0008         # Remote reactivation
  sigma    <- 1/0.5          # symptom development (6 mo)
  fast_p   <- 0.25           # Fraction fast progressing to active TB
  slow_p   <- 0.0008         # Remote reactivation
  sigma_p  <- 1/0.33         # symptom development (4 mo)
  imm      <- 0.5            # Infectiousness decline (partial immunity)
  I0       <- 1e-6
  P0       <- 250/1e5
  r_incar  <- 0.001          # rate of incarceration
  r_release<- 1/2.5          # rate of release (1/mean prison term)
  prob_inf_mix       <- 0.5  # probability of infection given short contact with external contactee
  case_detection     <- if (time > 950) 0.87 else 0#    # case detection rate
  case_detection_p   <- if (time > 950) 0.53 else 0#    # case detection rate
  tx_initiation_delay<- if (time > 950) 1/0.21 else 0#  # treatment initiation delay
  tx_iniation_delay_p<- if (time > 950) 1/0.27 else 0#  # treatment initiation delay
  tx_completion      <- if (time > 950) 0.8 else 0#     # fraction of people completing treatment 
  tx_duration        <- if (time > 950) 1/0.504 else 0# # duration treatment 184 days
  tx_forgiveness     <- if (time > 950) 0.5 else 0#     # fraction recovering although not completing treatment
  screen_delay       <- if (time > 950) 1/0.019 else 0# # delay of starting treatment after screening 7 days
  
  
  # External inputs ---------------------------------------------------------
  beta          <-parameter(5)
  beta_p        <-parameter(5)
  mix_p         <- parameter(0.008) # fraction of contacts from prison  (1 to 5% from Liu2024Lancet)
  screen_exit_a   <-parameter(0)
  screen_entry_a  <-parameter(0)
  screen_yearly_a <-parameter(0)
  screen_exit_s   <-parameter(0)
  screen_entry_s  <-parameter(0)
  screen_yearly_s <-parameter(0)
  
  # Initial conditions ------------------------------------------------
  
  
  initial(U) <- 1-I0-P0
  initial(L) <- 0
  initial(Ia) <- 0
  initial(Is) <- I0
  initial(Tr) <- 0
  initial(R) <- 0
  
  initial(U_p) <- P0-(I0*P0)
  initial(L_p) <- 0
  initial(Ia_p) <- 0
  initial(Is_p) <- I0*P0
  initial(Tr_p) <- 0
  initial(R_p) <- 0
  
  initial(incidence) <- 0 
  initial(incidence_p) <- 0 
  initial(new_imprisoned) <- 0
  initial(new_released)   <- 0
  initial(new_detected)   <- 0
  
  
})




# Process model output function -------------------------------------------


update_intervention<-function(mod,start_t,state0,run_t,param,label){
  
  # Set time   
  dust_system_set_time(mod, start_t)
  
  # Set model state  
  dust_system_set_state(mod,state0)
  
  # Set parameters
  dust_system_update_pars(sys, pars=param)
  
  # run simulation
  y0<- dust_system_simulate(mod, seq(start_t,start_t+run_t))
  
  # unpack results
  y0 <- dust_unpack_state(sys, y0)
  
  
  
  
  # Create output dataframe
  totalpop               = y0$Is + y0$Ia + y0$U + y0$L + y0$R + y0$Is_p + y0$Ia_p + y0$U_p + y0$L_p + y0$R_p + y0$Tr + y0$Tr_p
  prison_pop             = y0$Is_p + y0$Ia_p + y0$U_p + y0$L_p + y0$R_p + y0$Tr_p
  
  out<-data.frame(
    years                  = seq(0,run_t),  
    totalpop               = totalpop,
    prison_pop             = prison_pop,
    incidence_community    = y0$incidence,
    incidence_prison       = y0$incidence_p,
    n_annual_screened      = 6.8e6*(c(0,diff(y0$new_released+y0$new_imprisoned)) + prison_pop),
    n_cummulative_screened = cumsum(6.8e6*(c(0,diff(y0$new_released+y0$new_imprisoned)) + prison_pop)),
    n_annual_detected      = 6.8e6*(c(0,diff(y0$new_detected))),
    n_cumulative_detected  = cumsum(6.8e6*(c(0,diff(y0$new_detected)))),
    n_cummulative_tb_cases = 6.8e6*(y0$incidence + y0$incidence_p)
  )
  
  out$scenario<-label
  
  return(out)
  
}

