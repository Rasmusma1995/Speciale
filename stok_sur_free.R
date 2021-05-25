setwd("C:/Users/Rasmu/Desktop/speciale")
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(deSolve)
library(grid)
library(latex2exp)
library(furrr)
library(gridExtra)
library(ggpubr)


set.seed(100)

lf <- read_excel("C:/Users/Rasmu/Desktop/speciale/Benchmark_levetidsforbedringer2019 xls.xls")
obsdead <- read_excel("C:/Users/Rasmu/Desktop/speciale/Benchmark_doedelighed2019 xls.xls") 

alder <- 40
number_of_paths <- 1000
endtime <- 120

increment_size <- 0.1

obsdead_const <- obsdead[sum(floor(alder)+1),3] %>% pull
lf_constant <-  lf[sum(floor(alder)+1),3] %>% pull

#sim params 

b2 <-0.02
beta2 <- 0.02
sigma2 <- 0.15
sigma_af <- 0.10
X2_0 <- 1

#Technical intensities a x year old
tek_mu_01 <- function(t) (0.0004 + 10^(4.54+0.06*(t+alder)-10))*as.numeric(t+alder <= 65)
tek_mu_10 <- function(t) 2.0058*exp(-0.117*(t+alder))*as.numeric(t+alder <= 65)
tek_mu_02 <- function(t) 0.0005 + 10^(5.88 + 0.038*(t+alder)- 10)
tek_mu_12 <- function(t) tek_mu_02(t)*(1+ as.numeric(t+alder <= 65))

#markets intensties for a x year old

mu_01 <- function(t) as.numeric(t+alder <= 65)*(10^(5.662015+0.033462*(alder+t)-10))
mu_10 <- function(t) 4.0016*exp(-0.117*(alder+t))
mu_12 <- function(t) 0.010339+10^(5.070927+0.054049*(t+alder)-10)
mu_02 <- function(t) obsdead_const*(1-lf_constant)^t

mu_as_be <- function(t) 0.06-0.002*t
mu_af_be <- function(t) 0.05

#Payments 

t.r <- 0.01

b_i <- 100000
b0.65 <- 100000


# Utils -------------------------------------------------------------------

integrate_vec <- function(vec, mu_base, stepsize=increment_size, return_vec = T){
  
  time_grid <- seq(0, 65-alder,stepsize)
  tmp <- sapply(1:sum(length(vec)-1), function(i) vec[i]-mu_base(time_grid[i])+vec[i+1]-mu_base(time_grid[i+1]))
  
  if (return_vec) 
    return(cumsum(c(0,tmp/2*stepsize)))
  
  return(sum(tmp/2*stepsize))
}


insert_time_first_col <- 
  function(data, stepsize=increment_size,TotalTime=65-alder, opret_col = F ){
    
    time_grid <-  seq(0,TotalTime,stepsize)  
    
    if(opret_col){
      
      return(cbind("Time" = time_grid, data) %>% as.data.frame())
      
    }
    data[,1] <- time_grid
    colnames(data)[1] <- "Time"
    
    return(data)
  }

scaling_int_premium_to_surrender <- function(vec_sur, vec_free, mu_base, 
                                             stepsize=increment_size, returnWTime = F,rho = function(t) 1){
  
  time_grid <- seq(0, 65-alder,stepsize)
  
  
  rho_vec <- rho(time_grid)
  sur_part <- exp(-integrate_vec(vec_sur,mu_base))
  exp_free <- exp(-integrate_vec(vec_free,function(t) 0))
  
  free_policy_part <- sapply(1:sum(length(vec_free)-1), function(i) exp_free[i]*vec_free[i]*rho_vec[i]
                             +exp_free[i+1]*vec_free[i+1]*rho_vec[i+1])
  
  free_policy_part2 <- cumsum(c(0,free_policy_part/2*stepsize))
  
  if (returnWTime)
    return(data.frame(time = time_grid, val = sur_part*free_policy_part2))
  
  return(sur_part*free_policy_part2)
}

# Simulated the stochastic surrender & freepolicy --------------------------------------

stoc_surrender_free_sim <- function(n, stepsize, TotalTime){
  
  length <- TotalTime/stepsize
  
  s.Norm <- function(n) N[n-1]
  
  out <- matrix(NA, length+1 ,n+1) %>% as.data.frame()
  out2 <- matrix(NA, length+1 ,n+1) %>% as.data.frame()
  timegrid <- seq(0,TotalTime, stepsize)
  
  options(warn=-1)
  
  for (j in 1:n){
    sur_path <- X2_0
    free_path <- X2_0
    N <- rnorm(length,0, 1)
    for (i in 2:(length+1)){
      #surrender
      
      dsur <- (b2-beta2*sur_path)*stepsize+sigma2*sqrt(sur_path*stepsize)*s.Norm(i)
      sur_path <- sur_path + dsur
      out[i,j+1] <- mu_as_be(timegrid[i])*sur_path
      
      #free policy
      dfree <- (b2-beta2*free_path)*stepsize+sigma_af*sqrt(free_path*stepsize)*s.Norm(i)
      free_path <- free_path + dfree
      out2[i,j+1] <- mu_af_be(timegrid[i])*1/free_path
      
    }
    if (any(is.na(out[,j+1])) || any(is.na(out2[,j+1])))
      j <- j - 1
  }
  
  options(warn=0)
  #surrender
  
  out[1,] <- mu_as_be(0)
  out[,1] <- timegrid
  
  colnames(out)[1] <- "time"
  out$time <- as.character(out$time) %>% as.double()
  
  #free policy
  out2[1,] <- mu_af_be(0)
  out2[,1] <- timegrid
  
  colnames(out2)[1] <- "time"
  out2$time <- as.character(out2$time) %>% as.double()
  
  cat("NA detected in surrender sim: ", any(is.na(out)))
  cat("NA detected in free-policy sim: ", any(is.na(out2)))
  
  return(list(out=out, out2=out2))
  
}

# tekniske reserver -------------------------------------------------------

premium.cal <- function(t1, t2,finhed){
  
  odefunc <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dP00 <-  P01*tek_mu_10(t)- P00*(tek_mu_01(t) + tek_mu_02(t))
      dP01 <- P00*tek_mu_01(t)-P01*(tek_mu_12(t)+ tek_mu_10(t))
      dP11 <- P10*tek_mu_01(t)-P11*(tek_mu_10(t) + tek_mu_12(t))
      dP10 <- P11*tek_mu_10(t)-P10*(tek_mu_01(t) + tek_mu_02(t) )
      list(c(dP00,dP01, dP11, dP10))
    })
  }
  
  state      <- c(P00=1, P01=0, P11 = 1, P10=0)
  times      <- seq(t1, t2 , by = finhed)
  
  tek.probs <- ode(y = state, times = times, func = odefunc, parms = list(), method="rk4") %>% as.data.frame()
  
  l <- 1/finhed*(65-alder)+1
  
  premium <-tek.probs %>%
    mutate(numorator1 = ifelse(time<=65-alder, exp(-t.r*time)*P01*b_i, 0),
           numorator2 = ifelse(time>=65-alder, exp(-t.r*time)*(P01+P00)*b0.65, 0),
           denominator = ifelse(time<=65-alder, exp(-t.r*time)*P00, 0)) %>%
    dplyr::summarise(numorator1 = integrate_vec(.$numorator1[1:l], function(t) 0, return_vec = F),
                     numorator2 = integrate_vec(.$numorator2[-(1:(l-1))], function(t) 0, return_vec = F),
                     denominator= integrate_vec(.$denominator[1:l], function(t) 0, return_vec = F)) %>%
    mutate(p = (numorator1+numorator2)/denominator) %>%
    pull
  
  
  return(premium)
  
}

premium <- premium.cal(0, endtime-alder, 0.1);premium

tek_reserve_cal <- function(t1, t2,finhed, benefitsOnly = F){
  
  thiele <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dV0 <- t.r*V0 + premium*(!benefitsOnly)*(t<=65-alder)-b0.65*(t>=65-alder) - tek_mu_01(t)*(V1-V0)+tek_mu_02(t)*V0
      dV1 <- t.r*V1 - b_i -tek_mu_10(t)*(V0-V1)+tek_mu_12(t)*V1
      list(c(dV0,dV1))
    })
  }
  
  state      <- c(V0=0, V1=0)
  times      <- seq(t2, t1 , by = -finhed)
  
  out <- ode(y = state, times = times, func = thiele, parms = list(), method="rk4") %>% as.data.frame()
  
  return(out)
}


tek_reserve_to_rho <- tek_reserve_cal(0, endtime-alder, increment_size)
tek_reserve_to_rho_ben <- tek_reserve_cal(0, endtime-alder, increment_size, benefitsOnly=T)

rho_func <- approxfun(tek_reserve_to_rho[,1], tek_reserve_to_rho[,2]/tek_reserve_to_rho_ben[,2]) 

# solver for transition probabilities -------------------------------------

probs_market <- function(t1, t2,finhed, mu_base, mu_base_free){
  
  odefunc <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dP00 <-  P01*mu_10(t)- P00*(mu_base(t)+mu_01(t) + mu_02(t)+mu_base_free(t))
      dP01 <- P00*mu_01(t)-P01*(mu_12(t)+ mu_10(t))
      dP11 <- P10*mu_01(t)-P11*(mu_10(t) + mu_12(t))
      dP10 <- P11*mu_10(t)-P10*(mu_base(t)+mu_01(t) + mu_02(t) )
      dP04 <- P00*mu_base_free(t)*rho_func(t)+P05*mu_10(t) - P04*(mu_base(t)+mu_01(t) + mu_02(t))
      dP05 <- P04*mu_01(t) -P05*(mu_12(t)+ mu_10(t))
      dP14 <- P10*mu_base_free(t)*rho_func(t) + P15*mu_10(t)- P14*(mu_base(t)+mu_01(t) + mu_02(t))
      dP15 <- P14*mu_01(t)-P15*(mu_12(t)+ mu_10(t))
      list(c(dP00,dP01, dP11, dP10,dP04, dP05, dP14, dP15))
    })
  }
  
  state      <- c(P00=1, P01=0, P11 = 1, P10=0,P04=0, P05=0, P14 =0, P15=0)
  times      <- seq(t1, t2 , by = finhed)
  
  out <- ode(y = state, times = times, func = odefunc, parms = list(), method="rk4")
  
  return(out %>% as.data.frame())
  
}


# sim ---------------------------------------------------------------------

sim_data <- stoc_surrender_free_sim(number_of_paths, increment_size, 65-alder)

surrender_dat <- sim_data$out
free_policy_dat <- sim_data$out2

# mu_baseS = 0 -------------------------------------------------------------

transprops_mu_base_0 <- list(
  P00 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P01 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P11 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>%
    as.data.frame() %>%
    insert_time_first_col(),
  P10 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P04 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P05 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P14 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P15 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col()
)

#transition probabilites mu only be solved once. 

ode_mu_base_0 <- probs_market(0,65-alder, increment_size, function(t) 0, function(t) 0)

for (i in 2:ncol(surrender_dat)){
  transprops_mu_base_0$P00[,i] <- exp(-integrate_vec(surrender_dat[,i]+free_policy_dat[,i],mu_base = function(t) 0 ))*ode_mu_base_0[,"P00"]
  transprops_mu_base_0$P01[,i] <- exp(-integrate_vec(surrender_dat[,i]+free_policy_dat[,i],mu_base = function(t) 0))*ode_mu_base_0[,"P01"]
  transprops_mu_base_0$P11[,i] <- exp(-integrate_vec(surrender_dat[,i]+free_policy_dat[,i],mu_base = function(t) 0))*ode_mu_base_0[,"P11"]
  transprops_mu_base_0$P10[,i] <- exp(-integrate_vec(surrender_dat[,i]+free_policy_dat[,i],mu_base = function(t) 0))*ode_mu_base_0[,"P10"]
  transprops_mu_base_0$P04[,i] <- scaling_int_premium_to_surrender(surrender_dat[,i], free_policy_dat[,i],mu_base= function(t) 0, rho = rho_func)*ode_mu_base_0[,"P00"]
  transprops_mu_base_0$P05[,i] <- scaling_int_premium_to_surrender(surrender_dat[,i], free_policy_dat[,i],mu_base= function(t) 0, rho = rho_func)*ode_mu_base_0[,"P01"]
  transprops_mu_base_0$P14[,i] <- scaling_int_premium_to_surrender(surrender_dat[,i], free_policy_dat[,i],mu_base= function(t) 0, rho = rho_func)*ode_mu_base_0[,"P10"]
  transprops_mu_base_0$P15[,i] <- scaling_int_premium_to_surrender(surrender_dat[,i], free_policy_dat[,i],mu_base= function(t) 0, rho = rho_func)*ode_mu_base_0[,"P11"]
}

#Test that userdefined integration is correct

func_wrap_to_exclude <- function(time_to_test){

  test <- approxfun(surrender_dat[,1],surrender_dat[,2])
  test2 <- approxfun(free_policy_dat[,1],free_policy_dat[,2])
  test_free <- function(t){
    
    exp(-integrate(test2, lower=0, upper = t, subdivisions=300)$value)*test2(t)
  }
  out1 <-integrate(Vectorize(test_free), 0,time_to_test,subdivisions=300,stop.on.error = F)$value*exp(-integrate(test,lower = 0, upper = time_to_test, subdivisions=300)$value)
  out2 <- scaling_int_premium_to_surrender(surrender_dat[,2], free_policy_dat[,2], function(t) 0, returnWTime = T) %>%
    filter(time ==time_to_test)
  
  cat("R implementation gives :", out1, "\n")
  cat("Own implementation gives :", out2[,2], "\n")
}


# mu_baseS = E[mu_as] ------------------------------------------------------

mu_base_mean <- mu_as_be

transprops_mu_base_mean <-  list(
  P00 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P01 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P11 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>%
    as.data.frame() %>%
    insert_time_first_col(),
  P10 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P04 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P05 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P14 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P15 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col()
)

#transition probabilites mu only be solved once. 

ode_mu_base_mean <- probs_market(0,65-alder, increment_size, mu_base_mean, function(t) 0)

for (i in 2:ncol(surrender_dat)){
  transprops_mu_base_mean$P00[,i] <- exp(-integrate_vec(surrender_dat[,i]+free_policy_dat[,i],mu_base = mu_base_mean ))*ode_mu_base_mean[,"P00"]
  transprops_mu_base_mean$P01[,i] <- exp(-integrate_vec(surrender_dat[,i]+free_policy_dat[,i],mu_base = mu_base_mean))*ode_mu_base_mean[,"P01"]
  transprops_mu_base_mean$P11[,i] <- exp(-integrate_vec(surrender_dat[,i]+free_policy_dat[,i],mu_base = mu_base_mean))*ode_mu_base_mean[,"P11"]
  transprops_mu_base_mean$P10[,i] <- exp(-integrate_vec(surrender_dat[,i]+free_policy_dat[,i],mu_base = mu_base_mean))*ode_mu_base_mean[,"P10"]
  transprops_mu_base_mean$P04[,i] <- scaling_int_premium_to_surrender(surrender_dat[,i], free_policy_dat[,i],mu_base = mu_base_mean, rho = rho_func)*ode_mu_base_mean[,"P00"]
  transprops_mu_base_mean$P05[,i] <- scaling_int_premium_to_surrender(surrender_dat[,i], free_policy_dat[,i],mu_base = mu_base_mean, rho = rho_func)*ode_mu_base_mean[,"P01"]
  transprops_mu_base_mean$P14[,i] <- scaling_int_premium_to_surrender(surrender_dat[,i], free_policy_dat[,i],mu_base = mu_base_mean, rho = rho_func)*ode_mu_base_mean[,"P10"]
  transprops_mu_base_mean$P15[,i] <- scaling_int_premium_to_surrender(surrender_dat[,i], free_policy_dat[,i],mu_base = mu_base_mean, rho = rho_func)*ode_mu_base_mean[,"P11"]
}


# true model --------------------------------------------

transprops_mu_base_true <-  list(
  P00 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P01 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P11 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>%
    as.data.frame() %>%
    insert_time_first_col(),
  P10 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P04 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P05 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P14 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col(),
  P15 = matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>% 
    as.data.frame() %>%
    insert_time_first_col()
)

for (i in 2:ncol(surrender_dat)){
  mu_basetrue <-  approxfun(surrender_dat[,1], surrender_dat[,i]) 
  mu_basefree <- approxfun(free_policy_dat[,1], free_policy_dat[,i]) 
  ode_mu_base_true <- probs_market(0,65-alder, increment_size, mu_basetrue, mu_basefree)
  
  transprops_mu_base_true$P00[,i] <- ode_mu_base_true[,"P00"]
  transprops_mu_base_true$P01[,i] <- ode_mu_base_true[,"P01"]
  transprops_mu_base_true$P11[,i] <- ode_mu_base_true[,"P11"]
  transprops_mu_base_true$P10[,i] <- ode_mu_base_true[,"P10"]
  
  transprops_mu_base_true$P04[,i] <- ode_mu_base_true[,"P04"]
  transprops_mu_base_true$P05[,i] <- ode_mu_base_true[,"P05"]
  transprops_mu_base_true$P14[,i] <- ode_mu_base_true[,"P14"]
  transprops_mu_base_true$P15[,i] <- ode_mu_base_true[,"P15"]
}




# portfolio-wide means ----------------------------------------------------


v_stjerne_bar_circ <- function(sti = NULL, diff=F){
  
  
  tek_reserve <- tek_reserve_cal(0, endtime-alder, increment_size) %>%
    filter(time<=65-alder)
  tek_reserve_ben <- tek_reserve_cal(0, endtime-alder, increment_size) %>%
    filter(time<=65-alder)
  
  tmp_func <- function(ssh, sti){
    if (is.null(sti)){
      
      true <- rowMeans(transprops_mu_base_true[[ssh]][,-1]) 
      mean <- rowMeans(transprops_mu_base_mean[[ssh]][,-1]) 
      simpel <- rowMeans(transprops_mu_base_0[[ssh]][,-1]) 
      
    }
    else{
      true <- transprops_mu_base_true[[ssh]][,sum(sti+1)]
      mean <- transprops_mu_base_mean[[ssh]][,sum(sti+1)]
      simpel <- transprops_mu_base_0[[ssh]][,sum(sti+1)]
    }
    return(list(true=true, mean=mean, simpel=simpel))
  }
  
  P00 <- tmp_func("P00", sti)
  P01 <- tmp_func("P01", sti)
  P11 <- tmp_func("P11", sti)
  P10 <- tmp_func("P10", sti)
  P04 <- tmp_func("P04", sti)
  P05 <- tmp_func("P05", sti)
  P14 <- tmp_func("P14", sti)
  P15 <- tmp_func("P15", sti)
  
  V_z0_0_func <- function(char) P00[[char]]*rev(tek_reserve$V0)+P01[[char]]*rev(tek_reserve$V1) + P04[[char]]*rev(tek_reserve_ben$V0)+P05[[char]]*rev(tek_reserve_ben$V1)
  V_z0_1_func <- function(char) P11[[char]]*rev(tek_reserve$V1)+P10[[char]]*rev(tek_reserve$V0)+ P14[[char]]*rev(tek_reserve_ben$V0)+P15[[char]]*rev(tek_reserve_ben$V1)
  
  #z0 =0
  V_z0_0_true <- V_z0_0_func("true")
  V_z0_0_mean <- V_z0_0_func("mean")
  V_z0_0_simpel <- V_z0_0_func("simpel")
  
  #z0=1
  V_z0_1_true <- V_z0_1_func("true")
  V_z0_1_mean <- V_z0_1_func("mean")
  V_z0_1_simpel <- V_z0_1_func("simpel")
  
  if(diff){
    out <-  cbind(time = rev(tek_reserve$time),V_z0_0_true-V_z0_0_mean, 
                  V_z0_0_true-V_z0_0_simpel,
                  V_z0_1_true-V_z0_1_mean,V_z0_1_true-V_z0_1_simpel ) %>%
      as.data.frame()
    colnames(out) <- c("time", "V_z0_0_mean", "V_z0_0_simpel" 
                       ,"V_z0_1_mean", "V_z0_1_simpel")
    return(out)
  }
  
  out <- cbind(time = rev(tek_reserve$time),V_z0_0_true,V_z0_0_mean, V_z0_0_simpel,
               V_z0_1_true,V_z0_1_mean,V_z0_1_simpel ) %>%
    as.data.frame()
  
  return(out)
  
}


v_stjerne_bar_circ_plot_function <- 
  function(sti = NULL){
    
    time_grid <- seq(0, 65-alder, increment_size)
    
    tmp <- v_stjerne_bar_circ(sti, T)
    tmp2 <- v_stjerne_bar_circ(sti, F) %>%
      mutate(relv0_m = (V_z0_0_true-V_z0_0_mean)/V_z0_0_true*100,
             relv0_s = (V_z0_0_true-V_z0_0_simpel)/V_z0_0_true*100,
             relv1_m = (V_z0_1_true-V_z0_1_mean)/V_z0_1_true*100,
             relv1_s = (V_z0_1_true-V_z0_1_simpel)/V_z0_1_true*100) %>%
      select(time, relv0_m, relv0_s,relv1_m, relv1_s)
    
    p0 <- ggplot(data.frame(time_grid),aes(time_grid)) +
      geom_line(aes(y=tmp$V_z0_0_mean, color = "mean")) +
      geom_line(aes(y=tmp$V_z0_0_simpel, color = "simple")) +
      xlab("") + ylab(TeX("Expected difference")) + theme_bw() + ggtitle(TeX("\\{z_0 = 0\\}"))+
      scale_colour_manual(values = c('mean' = '#CC0000','simple' = 'black'),name = "", 
                          labels = expression("mean",0))
    
    p1 <- ggplot(data.frame(time_grid),aes(time_grid)) +
      geom_line(aes(y=tmp$V_z0_1_mean, color = "mean")) +
      geom_line(aes(y=tmp$V_z0_1_simpel, color = "simple")) +
      xlab("") + ylab(TeX("Expected difference")) + theme_bw() + ggtitle(TeX("\\{z_0 = 1\\}"))+
      scale_colour_manual(values = c('mean' = '#CC0000','simple' = 'black'),name = "", 
                          labels = expression("mean",0))
    
    p2 <- ggplot(data.frame(time_grid),aes(time_grid)) +
      geom_line(aes(y=tmp2$relv0_m, color = "mean")) +
      geom_line(aes(y=tmp2$relv0_s, color = "simple")) +
      xlab("") + ylab(TeX("Relative difference in %")) + theme_bw() + ggtitle(TeX("\\{z_0 = 0\\}"))+
      scale_colour_manual(values = c('mean' = '#CC0000','simple' = 'black'),name = "", 
                          labels = expression("mean",0))
    
    p3 <- ggplot(data.frame(time_grid),aes(time_grid)) +
      geom_line(aes(y=tmp2$relv1_m, color = "mean")) +
      geom_line(aes(y=tmp2$relv1_s, color = "simple")) +
      xlab("") + ylab(TeX("Relative difference in %")) + theme_bw() + ggtitle(TeX("\\{z_0 = 1\\}"))+
      scale_colour_manual(values = c('mean' = '#CC0000','simple' = 'black'),name = "", 
                          labels = expression("mean",0))
    
    
    
    p <- ggarrange(p0 +  theme(legend.position = "none"),
                   p2+ theme(legend.position = "none"),
                      p1+ theme(legend.position = "none"),
                      p3+ theme(legend.position = "none"),
                      ncol = 2,
                      nrow = 2,
                   common.legend = T,
                   legend = "top")
    
  }


vplot <- v_stjerne_bar_circ_plot_function()

annotate_figure(vplot,
                bottom = text_grob("Time", color = "black", vjust=-1, hjust=0))
