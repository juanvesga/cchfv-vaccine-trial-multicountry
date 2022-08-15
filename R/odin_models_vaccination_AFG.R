# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script contains the Odin compartmental models for CCHF in livestock (SIR)
# and among humans (SEIR). Also wrapper functions to format and package results
# to be used elsewhere


# Call Odin
library(odin)

## Deterministic SIR livestock model for CCHF transmission

sir_model <- odin::odin(
  {
  #ODEs
  deriv(S1) <- birth * (1 - prev) + Ri * lact_imm_loss + R1 * adult_imm_loss -
    S1 * foi - S1 * (mu[1] + ageing[1] + vacc_rate*vacc_eff)
  deriv(S2) <- S1 * ageing[1]  + R2 * adult_imm_loss - S2 * foi - S2 * (mu[2] + ageing[2]+ vacc_rate*vacc_eff)
  deriv(S3) <- S2 * ageing[2]  + R3 * adult_imm_loss - S3 * foi - S3 * (mu[3] + ageing[3]+ vacc_rate*vacc_eff)
  deriv(S4) <- S3 * ageing[3]  + R4 * adult_imm_loss - S4 * foi - S4 * (mu[4] + ageing[4]+ vacc_rate*vacc_eff)
  deriv(S5) <- S4 * ageing[4]  + R5 * adult_imm_loss - S5 * foi - S5 * (mu[5] + vacc_rate*vacc_eff) 
  
  deriv(V1) <- S1 * vacc_rate*vacc_eff - V1 * (mu[1] + ageing[1] + foi + time_protection)  
  deriv(V2) <- S2 * vacc_rate*vacc_eff + V1 * ageing[1] - V2 * (mu[2] + ageing[2] + foi + time_protection) 
  deriv(V3) <- S3 * vacc_rate*vacc_eff + V2 * ageing[2] - V3 * (mu[3] + ageing[3] + foi + time_protection)
  deriv(V4) <- S4 * vacc_rate*vacc_eff + V3 * ageing[3] - V4 * (mu[4] + ageing[4] + foi + time_protection)
  deriv(V5) <- S5 * vacc_rate*vacc_eff + V4 * ageing[4] - V5 * (mu[5] + foi + time_protection)
  
  deriv(Vp1) <- V1 * time_protection - Vp1 * (mu[1] + ageing[1])  
  deriv(Vp2) <- V2 * time_protection + Vp1 * ageing[1] - Vp2 * (mu[2] + ageing[2])
  deriv(Vp3) <- V3 * time_protection + Vp2 * ageing[2]- Vp3 * (mu[3] + ageing[3])
  deriv(Vp4) <- V4 * time_protection + Vp3 * ageing[3]- Vp4 * (mu[4] + ageing[4])
  deriv(Vp5) <- V5 * time_protection + Vp4 * ageing[4]- Vp5 * (mu[5]) 
  
  
  deriv(I1) <- (S1 + V1) * foi + seed[1] - I1 * (recovery + mu[1] + ageing[1])
  deriv(I2) <- (S2 + V2) * foi + seed[2] + I1 * ageing[1] - I2 * (recovery + mu[2] + ageing[2])
  deriv(I3) <- (S1 + V3) * foi + seed[3] + I2 * ageing[2] - I3 * (recovery + mu[3] + ageing[3])
  deriv(I4) <- (S4 + V4) * foi + seed[4] + I3 * ageing[3] - I4 * (recovery + mu[4] + ageing[4])
  deriv(I5) <- (S5 + V5) * foi + seed[5] + I4 * ageing[4] - I5 * (recovery + mu[5])
  
  deriv(Ri) <- birth * prev - Ri * (lact_imm_loss + mu[1])
  deriv(R1) <- I1 * recovery - R1 * (mu[1] + ageing[1] + adult_imm_loss)
  deriv(R2) <- I2 * recovery + R1 * ageing[1] - R2 * (mu[2] + ageing[2] + adult_imm_loss)
  deriv(R3) <- I3 * recovery + R2 * ageing[2] - R3 * (mu[3] + ageing[3] + adult_imm_loss)
  deriv(R4) <- I4 * recovery + R3 * ageing[3] - R4 * (mu[4] + ageing[4] + adult_imm_loss)
  deriv(R5) <- I5 * recovery + R4 * ageing[4] - R5 * (mu[5] + adult_imm_loss)
  
  output(vaccinated)<- (S1+S2+S3+S4+S5)*vacc_rate
  
  # initial conditions
  initial(S1) <- states0[1]
  initial(S2) <- states0[2]
  initial(S3) <- states0[3]
  initial(S4) <- states0[4]
  initial(S5) <- states0[5]
  initial(V1) <- states0[6]
  initial(V2) <- states0[7]
  initial(V3) <- states0[8]
  initial(V4) <- states0[9]
  initial(V5) <- states0[10]
  initial(Vp1)<- states0[11]
  initial(Vp2)<- states0[12]
  initial(Vp3)<- states0[13]
  initial(Vp4)<- states0[14]
  initial(Vp5)<- states0[15]
  initial(I1) <- states0[16]
  initial(I2) <- states0[17]
  initial(I3) <- states0[18]
  initial(I4) <- states0[19]
  initial(I5) <- states0[20]
  initial(Ri) <- states0[21]
  initial(R1) <- states0[22]
  initial(R2) <- states0[23]
  initial(R3) <- states0[24]
  initial(R4) <- states0[25]
  initial(R5) <- states0[26]
  
  # Necessary variables 
  N <- S1 + S2 + S3 + S4 + S5 +
    V1 + V2 + V3 + V4 + V5 + 
    Vp1 + Vp2 + Vp3 + Vp4 + Vp5 + 
    I1 + I2 + I3 + I4 + I5 + 
    Ri + R1 + R2 + R3 + R4 + R5 
  
  infectious <- I1 + I2 + I3 + I4 + I5
  
  birth <- (S1+I1+V1+Vp1+R1+Ri) * mu[1] + 
    (S2+V2+Vp2+I2+R2)    * mu[2] +
    (S3+V3+Vp3+I3+R3)    * mu[3] +
    (S4+V4+Vp4+I4+R4)    * mu[4] +
    (S5+V5+Vp5+I5+R5)    * mu[5]  
  
  prev <- (R1 + R2 + R3 + R4 + R5 + Ri)/N
  
  # transovarian seed every year
  
  seed[1]<- seed_t[1] * ((t%%12)==0)
  seed[2]<- seed_t[2] * ((t%%12)==0)
  seed[3]<- seed_t[3] * ((t%%12)==0)
  seed[4]<- seed_t[4] * ((t%%12)==0)
  seed[5]<- seed_t[5] * ((t%%12)==0)
  
  # Transition rates and user defined 
  envr_beta_vec <- interpolate(envr_beta_t, envr_beta_y, "linear")
  vacc_rate     <- interpolate(envr_beta_t, newvacc_y, "linear")
  foi <- 1- exp(-((envr_beta_vec/duration)*infectious/N))

  
  mu[1]            <- params[1]
  mu[2]            <- params[2]
  mu[3]            <- params[3]
  mu[4]            <- params[4]
  mu[5]            <- params[5]
  ageing[1]        <- params[6]
  ageing[2]        <- params[7]
  ageing[3]        <- params[8]
  ageing[4]        <- params[9]
  ageing[5]        <- params[10]
  lact_imm_loss   <- params[11] 
  adult_imm_loss  <- params[12]
  recovery        <- params[13]
  time_protection <- params[14]
  vacc_eff        <- params[15]
  duration        <-  1/recovery
  
  states0[]     <- user()
  envr_beta_t[] <- user()
  envr_beta_y[] <- user()
  params[]      <- user()
  seed_t[]      <- user()
  newvacc_y[]   <- user()
  
  # boilerplate from odin 
  dim(states0)    <-user()
  dim(envr_beta_t)<-user()
  dim(envr_beta_y)<-user()
  dim(params)     <-user()
  dim(seed_t)     <-user()
  dim(newvacc_y)  <-user()
  dim(seed)       <-5
  dim(mu)         <-5
  dim(ageing)     <-5
  
  
})



## Stochastic SEIR for spillover from livestock into humans

seir_stoc <- odin::odin(
  {
  ## Core equations for transitions between compartments:
  update(F_S) <- F_S - newInfections_F - newVacc_F_eff + newSusceptible_F + Birth_F - newdFS
  update(F_V) <- F_V + newVacc_F_eff - newVaccImm_F - newInfectionsV_F - newdFV
  update(F_Vp)<- F_Vp + newVaccImm_F - newdFVp
  update(F_E) <- F_E + newInfections_F + newInfectionsV_F - newInfectious_F  - newdFE
  update(F_I) <- F_I + newInfectious_F - newRecovered_F   - newdFI - newfatal_F
  update(F_R) <- F_R + newRecovered_F - newSusceptible_F  - newdFR
  
  # Infection in Others
  update(O_S) <- O_S - newInfections_O - newVacc_O_eff + newSusceptible_O + Birth_O - newdOS
  update(O_V) <- O_V + newVacc_O_eff - newVaccImm_O - newInfectionsV_O - newdOV
  update(O_Vp)<- O_Vp+ newVaccImm_O - newdOVp
  update(O_E) <- O_E + newInfections_O + newInfectionsV_O - newInfectious_O  - newdOE
  update(O_I) <- O_I + newInfectious_O - newRecovered_O   - newdOI - newfatal_O
  update(O_R) <- O_R + newRecovered_O - newSusceptible_O  - newdOR
  
  output(F_incidence) <- newInfectious_F
  output(O_incidence) <- newInfectious_O
  output(F_fatality)  <- newfatal_F 
  output(O_fatality)  <- newfatal_O 
  output(vaccinated)  <- newVacc_F + newVacc_O 
  
  
  # Initial conditions
  initial(F_S) <- states0[1]
  initial(F_V) <- states0[2]
  initial(F_Vp)<- states0[3]
  initial(F_E) <- states0[4]
  initial(F_I) <- states0[5]
  initial(F_R) <- states0[6]
  initial(O_S) <- states0[7]
  initial(O_V) <- states0[8]
  initial(O_Vp)<- states0[9]
  initial(O_E) <- states0[10]
  initial(O_I) <- states0[11]
  initial(O_R) <- states0[12]
  
  ################### HUMANS
  ProbaInf_F <- 1 - exp(-beta_farmer*risk_livestock)
  ProbaInf_O <- 1 - exp(-beta_other*risk_livestock)
  p_VtoVp    <- 1 - exp(-time_protection)
  p_EtoI     <- 1 - exp(-time_to_infous)
  p_ItoR     <- 1 - exp(-time_immune_human)
  p_RtoS     <- 1 - exp(-time_susceptible_human)
  p_bd       <- 1 - exp(-b_d)
  p[1] <- 1 - CFR*clin_frac
  p[2] <- CFR*clin_frac
  
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  
  # farmers S--> E: Prob of infection
  nS_outF         <- rpois(lambda=F_S*(ProbaInf_F+p_bd))
  newInfections_F <- round(rbinom(nS_outF,ProbaInf_F/(ProbaInf_F+p_bd)))
  newVacc_F       <- round(F_S*vaccRate_F)
  newVacc_F_eff   <- round(F_S*vaccRate_F*vacc_eff)
  
  # V to E
  newInfectionsV_F<- rpois(lambda=F_V*(ProbaInf_F))
  nS_outVf       <- round(rbinom((F_V-newInfectionsV_F),(p_VtoVp+p_bd)))
  # V to Vp
  newVaccImm_F    <- round(rbinom(nS_outVf,p_VtoVp/(p_VtoVp+p_bd)))
  
  # E-->I
  newInfectious_F <- round(F_E*p_EtoI)
  # I --> R
  noutI_F_all     <- round(rbinom(F_I,(p_ItoR+p_bd)))
  noutI_F         <- round(rbinom(noutI_F_all,p_ItoR/(p_ItoR+p_bd)))
  outcome_F[]     <- rmultinom(noutI_F,p)
  newRecovered_F  <- outcome_F[1]
  newfatal_F      <- outcome_F[2]
  
  # R --> S
  newSusceptible_F <- round(F_R*p_RtoS)
  
  # OTHERS S--> E: Prob of infection
  nS_outO <-  rpois(lambda=O_S*(ProbaInf_O+p_bd))
  newInfections_O <- round(rbinom(nS_outO,ProbaInf_O/(ProbaInf_O+p_bd)))
  newVacc_O       <- round(O_S*vaccRate_O)
  newVacc_O_eff   <- round(O_S*vaccRate_O*vacc_eff)
  
  # V E
  newInfectionsV_O<- rpois(lambda=O_V*(ProbaInf_O))
  nS_outVo       <- round(rbinom((O_V-newInfectionsV_O),(p_VtoVp+p_bd)))
  # V to Vp
  newVaccImm_O    <- round(rbinom(nS_outVo,p_VtoVp/(p_VtoVp+p_bd)))
  # E-->I
  newInfectious_O <- round(O_E*p_EtoI)
  # I --> R
  noutI_O_all     <- round(rbinom(O_I,(p_ItoR+p_bd)))
  noutI_O         <- round(rbinom(noutI_O_all,p_ItoR/(p_ItoR+p_bd)))
  outcome_O[]     <- rmultinom(noutI_O,p)
  newRecovered_O  <- outcome_O[1]
  newfatal_O      <- outcome_O[2]
  # R --> S
  newSusceptible_O <- round(O_R*p_RtoS)
  
  newdFS <- nS_outF - newInfections_F 
  newdFV <- round(F_V*p_bd)
  newdFVp<- nS_outVf-newVaccImm_F 
  newdFE <- round(F_E*p_bd)
  newdFI <- noutI_F_all - noutI_F 
  newdFR <- round(F_R*p_bd)
  newdOS <- nS_outO - newInfections_O 
  newdOV <- round(O_V*p_bd)
  newdOVp<- nS_outVo-newVaccImm_O 
  newdOE <- round(O_E*p_bd)
  newdOI <- noutI_O_all - noutI_O 
  newdOR <- round(O_R*p_bd)
  # 
  Birth_F <- newdFS+newdFE+newdFI+newdFR+newdFV+newdFVp
  Birth_O <- newdOS+newdOE+newdOI+newdOR+newdOV+newdOVp
  
  # User defined input
  
  # Transition rates and user defined 
  risk_livestock <- interpolate(risk_livestock_t, risk_livestock_y, "linear")
  vaccRate_F <- interpolate(risk_livestock_t, n_vaccination_f_y, "linear")
  vaccRate_O <- interpolate(risk_livestock_t, n_vaccination_o_y, "linear")
  
  b_d              <- params[1]
  beta_farmer      <- params[2]
  beta_other       <- params[3]
  time_to_infous   <- params[4]
  time_immune_human<- params[5]
  time_susceptible_human <- params[6]
  time_protection  <- params[7]
  CFR              <- params[8] 
  clin_frac        <- params[9]
  vacc_eff         <- params[10]
  
  states0[]            <- user()
  risk_livestock_t[]   <- user()
  risk_livestock_y[]   <- user()
  n_vaccination_f_y[]  <- user()
  n_vaccination_o_y[]  <- user()
  params[]             <- user()
  
  # boilerplate from odin 
  dim(states0)           <-user()
  dim(risk_livestock_t)  <-user()
  dim(risk_livestock_y)  <-user()
  dim(n_vaccination_f_y) <-user()
  dim(n_vaccination_o_y) <-user()
  dim(params)            <-user()
  dim(p) <- 2
  dim(outcome_F) <- 2
  dim(outcome_O) <- 2
  
  
}, verbose = FALSE)

# Wrapping function to get parameters , orgaanize input, create model and run SIR

run_sir <- function(params,theta)
  {
  
  ## Livestock SIR initial conditions
  init_prev<-sum(params$imm_t0*params$pop_st*params$N_liv)/params$N_liv
  prev_vec<-rep(0,5) 
  prev_vec[1]<-init_prev
  
  init.sir <- c(
    L_S = (1-prev_vec)*(1-params$imm_t0)*params$pop_st*params$N_liv - 
      params$pop_st,
    L_V = c(0,0,0,0,0),
    L_Vp= c(0,0,0,0,0),
    L_I = params$pop_st,
    L_Ri= prev_vec[1]*(1-params$imm_t0[1])*params$pop_st[1]*params$N_liv ,
    L_R = params$imm_t0*params$pop_st*params$N_liv
  )
  

  # Calculating R value for livestock environmental dependent
  envr_data<-data.frame(
    times=seq(1,params$nt), 
    envr_factor=environment_factor)
  
  foi_envr_factor_df<-envr_data%>%
    group_by(times)%>%
    mutate(envr_R_fac= envrdriver_foi_func(envr_factor,params$theta[["A"]]))%>%
    select(times,envr_R_fac)
  
  envr_beta_t <- foi_envr_factor_df$times
  envr_beta_y <- foi_envr_factor_df$envr_R_fac
  
  sir.params      <- c(
    params$deathd,
    params$Ageing,
    params$time_passimm_loss_livestock,
    params$time_susceptible_livestock,
    params$recover,
    params$time_protection,
    params$vacc_eff
  )
  
  # Vaccinated numbers
  doses <- round(params$p_vacc_livestock*100)

  VacPerMonth<- rmultinom(n=1,size=doses,prob=rep(1,params$stop_vax-params$start_vax))
  
  n_vaccination_y <- seq(1,params$nt)*0
  vac_cov<-VacPerMonth/100
  campaign_times<- params$start_vax:(params$stop_vax-1)
  n_vaccination_y[campaign_times]<-vac_cov
  
  if(params$liv_vacc_yearly==TRUE){
    yrs<-floor((length(n_vaccination_y)-params$stop_vax)/12)
    
  
    for (jj in 1:length(vac_cov)){
      times_yearly<-seq(from=campaign_times[jj]+12, to=campaign_times[jj]+12*yrs,12)  
      n_vaccination_y[times_yearly]<-vac_cov[jj]
      
    }
  }

  

  # Generate a model instance with the input
  sir.instance<-sir_model$new(states0=init.sir, 
                              envr_beta_t=envr_beta_t, 
                              envr_beta_y=envr_beta_y,
                              params=sir.params,
                              seed_t=params$pop_st,
                              newvacc_y=n_vaccination_y)
  
  
  # Run the model 
  t <- seq(1,params$nt)
  out.sir <- sir.instance$run(t)
  
  return(out.sir)
  
}  


# Call model function and organize/format model outputs SIR

get_output_sir <- function(params,theta)
  {
  nruns<- nrow(theta)
  mo_length<-params$nt
  yr_length<-length(unique(as.Date(temp_month$month,format="%Y")))
  
  # Allocate memory
  l_prev_age <- matrix(NA, nrow=nruns, ncol=5)
  l_prev_age_peak <- matrix(NA, nrow=nruns, ncol=5)
  l_prev_age_low <- matrix(NA, nrow=nruns, ncol=5)
  l_prev_all <- matrix(NA, nrow=nruns, ncol=1)
  l_prev_all_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  n_vaccinated  <- matrix(NA, nrow=nruns, ncol=mo_length)
  l_prev_0_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  l_prev_1_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  l_prev_2_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  l_prev_3_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  l_prev_4_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  livestock_pop <- matrix(NA, nrow=nruns, ncol=mo_length)
  sus_l_frac <- matrix(NA, nrow=nruns, ncol=mo_length)
  l_prev_inf_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  
  
  for (jj in 1:nruns){
    
    
    # Pass calibrated parameters 
    params$theta <- c(A=theta$A[jj], 
                      D_lact_liv=theta$D_lact_liv[jj],
                      D_imm_liv=theta$D_imm_liv[jj])
    
    params$time_passimm_loss_livestock<- 1/params$theta[['D_lact_liv']] 
    params$time_susceptible_livestock<- 1/params$theta[['D_imm_liv']] 
    
    # Call model function 
    
    out<-as.data.frame(run_sir(params))
    
    #------------------------------------ 
    # Process and format model output 
    #------------------------------------ 
    
    
    # Livestock prevalence by age 
    L_1<-out$S1+out$V1+out$Vp1+out$I1+out$R1+out$Ri
    L_2<-out$S2+out$V2+out$Vp2+out$I2+out$R2
    L_3<-out$S3+out$V3+out$Vp3+out$I3+out$R3
    L_4<-out$S4+out$V4+out$Vp4+out$I4+out$R4
    L_5<-out$S5+out$V5+out$Vp5+out$I5+out$R5
    
    
    l_prev_0_long[jj,]<-c(out$Ri+out$R1)/L_1
    l_prev_1_long[jj,]<-out$R2/L_2
    l_prev_2_long[jj,]<-out$R3/L_3
    l_prev_3_long[jj,]<-out$R4/L_4
    l_prev_4_long[jj,]<-out$R5/L_5
    
    simu_all_N_long<-c(L_1+L_2+L_3+L_4+L_5)
    livestock_pop[jj,]<- simu_all_N_long 
    
    simu_age_R<-c(
      mean(out$R1[(field_work_start:field_work_end)]+out$Ri[(field_work_start:field_work_end)]),
      mean(out$R2[(field_work_start:field_work_end)]),
      mean(out$R3[(field_work_start:field_work_end)]),
      mean(out$R4[(field_work_start:field_work_end)]),
      mean(out$R5[(field_work_start:field_work_end)]))
    
    simu_age_N<-c(
      mean(L_1[(field_work_start:field_work_end)]),
      mean(L_2[(field_work_start:field_work_end)]),
      mean(L_3[(field_work_start:field_work_end)]),
      mean(L_4[(field_work_start:field_work_end)]),
      mean(L_5[(field_work_start:field_work_end)]))
    
    simu_age_R_peak<-c(
      mean(out$R1[(field_work_start+1)]+out$Ri[(field_work_start+1)]),
      mean(out$R2[(field_work_start+1)]),
      mean(out$R3[(field_work_start+1)]),
      mean(out$R4[(field_work_start+1)]),
      mean(out$R5[(field_work_start+1)]))
    
    simu_age_N_peak<-c(
      mean(L_1[(field_work_start+1)]),
      mean(L_2[(field_work_start+1)]),
      mean(L_3[(field_work_start+1)]),
      mean(L_4[(field_work_start+1)]),
      mean(L_5[(field_work_start+1)]))
    
    simu_age_R_low<-c(
      mean(out$R1[(field_work_start-1)]+out$Ri[(field_work_start-1)]),
      mean(out$R2[(field_work_start-1)]),
      mean(out$R3[(field_work_start-1)]),
      mean(out$R4[(field_work_start-1)]),
      mean(out$R5[(field_work_start-1)]))
    
    simu_age_N_low<-c(
      mean(L_1[(field_work_start-1)]),
      mean(L_2[(field_work_start-1)]),
      mean(L_3[(field_work_start-1)]),
      mean(L_4[(field_work_start-1)]),
      mean(L_5[(field_work_start-1)]))
    
    
    l_prev_age[jj,] <- simu_age_R / simu_age_N
    l_prev_age_peak[jj,] <- simu_age_R_peak / simu_age_N_peak
    l_prev_age_low[jj,] <- simu_age_R_low / simu_age_N_low
    
    simu_all_R<-c(out$Ri+out$R1+out$R2+out$R3+out$R4+out$R5)
    l_prev_all_long[jj,]<-simu_all_R/simu_all_N_long
    
    simu_all_R<-mean(simu_all_R[(field_work_start:field_work_end)])
    simu_all_N<-mean(simu_all_N_long[(field_work_start:field_work_end)])
    l_prev_all[jj]<-simu_all_R/simu_all_N
    l_prev_0_long[jj,]<-c(out$Ri+out$R1)/L_1
    l_prev_1_long[jj,]<-out$R2/L_2
    l_prev_2_long[jj,]<-out$R3/L_3
    l_prev_3_long[jj,]<-out$R4/L_4
    l_prev_4_long[jj,]<-out$R5/L_5
    
    
    
    # susceptible fraction of livestock
    sus_l_frac[jj,]<-c(out$S1+out$S2+out$S3+out$S4+
                         out$S5)/simu_all_N_long
    
    l_prev_inf_long[jj,]<-c(out$I1+out$I2+out$I3+out$I4+
                              out$I5)/simu_all_N_long
    
    n_vaccinated[jj,] <-out$vaccinated
    
  }
  
  
  
  sim<-list(
    
    l_prev_age = l_prev_age,
    l_prev_age_peak = l_prev_age_peak,
    l_prev_age_low = l_prev_age_low,
    l_prev_all = l_prev_all,
    l_prev_all_long = l_prev_all_long,
    l_prev_0_long=l_prev_0_long,
    l_prev_1_long=l_prev_1_long,
    l_prev_2_long=l_prev_2_long,
    l_prev_3_long=l_prev_3_long,
    l_prev_4_long=l_prev_4_long,
    sus_l_frac = sus_l_frac,
    prev_inf_liv=l_prev_inf_long,
    livestock_pop=livestock_pop,
    l_vaccinated=n_vaccinated
    
  )
  
  
  return(sim)
  
  
}

# Wrapping function to get parameters , orgaanize input, create model and run SEIR

run_seir <- function(init.seir,params,times,n.replicates=1)
  {
  
  risk_livestock_t <- seq(1,params$nt)
  risk_livestock_y <- params$prev_inf_liv
  
  
  seir.params<-c(
    params$b_d,  
    params$theta[["F_risk"]],
    params$theta[["F_risk"]]*params$theta[["O_factor"]],
    params$time_to_infous,
    params$time_immune_human,
    params$time_susceptible_human,
    params$time_protection,
    params$CFR,
    params$theta[["clin_frac"]],
    params$vacc_eff
    
  )
  
  # Vaccination vector 
  
  doses_f <- round(params$pop_coverage_f*100)
  doses_o <- round(params$pop_coverage_o*100)
  
  VacPerMonth_f <- rmultinom(n=1,size=doses_f,prob=rep(1,params$stop_vax-params$start_vax))
  
  
  VacPerMonth_o <- rmultinom(n=1,size=doses_o,prob=rep(1,params$stop_vax-params$start_vax))
  
  
  n_vaccination_f_y <- seq(1,params$nt)*0
  n_vaccination_o_y <- seq(1,params$nt)*0
  n_vaccination_f_y[params$start_vax:(params$stop_vax-1)]<-VacPerMonth_f/100
  n_vaccination_o_y[params$start_vax:(params$stop_vax-1)]<-VacPerMonth_o/100  
  
  
  # Generate a stoch model instance with the input
  seir.instance<-seir_stoc$new(states0=init.seir, 
                               risk_livestock_t=risk_livestock_t, 
                               risk_livestock_y=risk_livestock_y,
                               n_vaccination_f_y=n_vaccination_f_y,
                               n_vaccination_o_y=n_vaccination_o_y,
                               params=seir.params)
  
  
  # Run the model 
  t <- times
  
  
  
  
  if (n.replicates>1){
    
    out.seir<-(replicate(100, seir.instance$run(t)[, -1]))

    # Find mean of 100 replicates of stoch model  
    out<-rowMeans(out.seir,na.rm = TRUE, dims=2)
    
    
  } else{
    out.seir <- seir.instance$run(t)
    out<-out.seir
  }
  
  
  return(out)
  
}  



##########
## End of Code
##########
