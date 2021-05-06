setwd("C:/Users/Rasmu/Desktop/speciale")
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(deSolve)
library(furrr)
library(gridExtra)

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




# Simulated the stochastic surrender --------------------------------------

stoc_surrender_sim <- function(n, stepsize, TotalTime){
  
  length <- TotalTime/stepsize
  
  s.Norm <- function(n) N[n-1]
  
  out <- matrix(NA, length+1 ,n+1) %>% as.data.frame()
  timegrid <- seq(0,TotalTime, stepsize)
 
  options(warn=-1)
  
  for (j in 1:n){
    X2_path <- X2_0
    N <- rnorm(length,0, 1)
    for (i in 2:(length+1)){
      dX_2 <- (b2-beta2*X2_path)*stepsize+sigma2*sqrt(X2_path*stepsize)*s.Norm(i)
      X2_path <- X2_path + dX_2
      out[i,j+1] <- mu_as_be(timegrid[i])*X2_path
    }
    if (any(is.na(out[,j+1])))
      j <- j - 1
  }
  
  options(warn=0)
  
  out[1,] <- mu_as_be(0)
  out[,1] <- timegrid
  
  colnames(out)[1] <- "time"
  out$time <- as.character(out$time) %>% as.double()
  
  cat("NA detected in surrender sim: ", any(is.na(out)))
  
  return(out)
  
}


# solver for transition probabilities -------------------------------------

probs_market <- function(t1, t2,finhed, mu_base){
  
  odefunc <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dP00 <-  P01*mu_10(t)- P00*(mu_base(t)+mu_01(t) + mu_02(t))
      dP01 <- P00*mu_01(t)-P01*(mu_12(t)+ mu_10(t))
      dP11 <- P10*mu_01(t)-P11*(mu_10(t) + mu_12(t))
      dP10 <- P11*mu_10(t)-P10*(mu_base(t)+mu_01(t) + mu_02(t) )
      list(c(dP00,dP01, dP11, dP10))
    })
  }
  
  state      <- c(P00=1, P01=0, P11 = 1, P10=0)
  times      <- seq(t1, t2 , by = finhed)
  
  out <- ode(y = state, times = times, func = odefunc, parms = list(), method="rk4")
  
  return(out %>% as.data.frame())
  
}


# Sim surrender -----------------------------------------------------------

surrender_dat <- stoc_surrender_sim(number_of_paths, increment_size, 65-alder)


# mu_base = 0 -------------------------------------------------------------

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
    insert_time_first_col()
)

#transition probabilites mu only be solved once. 

ode_mu_base_0 <- probs_market(0,65-alder, increment_size, function(t) 0)

for (i in 2:ncol(surrender_dat)){
  transprops_mu_base_0$P00[,i] <- exp(-integrate_vec(surrender_dat[,i],mu_base = function(t) 0 ))*ode_mu_base_0[,"P00"]
  transprops_mu_base_0$P01[,i] <- exp(-integrate_vec(surrender_dat[,i],mu_base = function(t) 0))*ode_mu_base_0[,"P01"]
  transprops_mu_base_0$P11[,i] <- exp(-integrate_vec(surrender_dat[,i],mu_base = function(t) 0))*ode_mu_base_0[,"P11"]
  transprops_mu_base_0$P10[,i] <- exp(-integrate_vec(surrender_dat[,i],mu_base = function(t) 0))*ode_mu_base_0[,"P10"]
}

# mu_base = E[mu_as] ------------------------------------------------------

mu_base_mean <- mu_as_be

transprops_mu_base_mean <- list(
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
    insert_time_first_col()
)

#transition probabilites mu only be solved once. 

ode_mu_base_mean <- probs_market(0,65-alder, increment_size, mu_base_mean)

for (i in 2:ncol(surrender_dat)){
  transprops_mu_base_mean$P00[,i] <- exp(-integrate_vec(surrender_dat[,i], mu_base = mu_base_mean))*ode_mu_base_mean[,"P00"]
  transprops_mu_base_mean$P01[,i] <- exp(-integrate_vec(surrender_dat[,i],  mu_base = mu_base_mean))*ode_mu_base_mean[,"P01"]
  transprops_mu_base_mean$P11[,i] <- exp(-integrate_vec(surrender_dat[,i],  mu_base = mu_base_mean))*ode_mu_base_mean[,"P11"]
  transprops_mu_base_mean$P10[,i] <- exp(-integrate_vec(surrender_dat[,i],  mu_base = mu_base_mean))*ode_mu_base_mean[,"P10"]
}


# mu_base = forward -------------------------------------------------------

#Kalibrating the forward curve

mu_base_forward_tmp1 <- matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>%
  as.data.frame() %>%
  insert_time_first_col()

mu_base_forward_tmp2 <- matrix(NA, nrow=nrow(surrender_dat), ncol = ncol(surrender_dat)) %>%
  as.data.frame() %>%
  insert_time_first_col()



for (i in 2:ncol(surrender_dat)){
  mu_base_forward_tmp1[,i] <- exp(-integrate_vec(surrender_dat[,i], function(t) 0))
  mu_base_forward_tmp2[,i] <- surrender_dat[,i]*mu_base_forward_tmp1[,i]
}

mu_base_forward <- approxfun(mu_base_forward_tmp1[,"Time"], rowMeans(mu_base_forward_tmp2[,-1])/rowMeans(mu_base_forward_tmp1[,-1]))

# odeSolveAffine1 <- function(t1, t2,finhed){
# 
#   odefunc <- function(t, state, parameters) {
#     with(as.list(c(state, parameters)), {
#       dphi <- -b2*psi
#       dpsi <- -1/2*sigma2^2*psi^2+beta2*psi+mu_as_be(t)
#       list(c(dphi,dpsi))
#     })
#   }
# 
#   state      <- c(phi = 0, psi=0)
#   times      <- seq(t2, t1 , by = -finhed)
# 
#   out <- ode(y = state, times = times, func = odefunc, parms = list(), method="euler")%>% as.data.frame()
# 
# 
#   return(out)
# 
# }
# 
# odeSolveAffine2 <- function(t1,t2,finhed){
# 
#   solve_tmp <- odeSolveAffine1(t1,t2, finhed)
#   phi <- approxfun(solve_tmp[,1], solve_tmp[,2], method = "constant")
#   psi <- approxfun(solve_tmp[,1], solve_tmp[,3], method = "constant")
# 
# 
#   odefunc <- function(t, state, parameters) {
#     with(as.list(c(state, parameters)), {
#       dA <- -b2*B
#       dB <- -psi(t)*sigma2^2*B+beta2*B
#       list(c(dA,dB))
#     })
#   }
# 
#   state      <- c(A =0, B=mu_as_be(t2))
#   times      <- seq(t2, t1 , by = -finhed)
# 
#   out <- ode(y = state, times = times, func = odefunc, parms = list(), method="rk4")%>% as.data.frame()
# 
# }
# 
# forward_cal <- function(s, finhed, start){
# 
#   grid <- seq(finhed,s,finhed)
#   out <- c(mu_as_be(0))
#   for (i in grid){
#     tmp <- odeSolveAffine2(0, i,finhed) %>% tail(1)
#     out <- append(out, tmp$A+tmp$B*start)
#   }
# 
#   return(data.frame(time = c(0,grid), forward = out))
# }
# 
# forward_rate_sol <- forward_cal(25, 0.1, 1)
# mu_base_forward <- approxfun(forward_rate_sol[,"time"],forward_rate_sol[,"forward"])


#Transprop 

transprops_mu_base_forward <- list(
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
    insert_time_first_col()
)

#transition probabilites mu only be solved once. 

ode_mu_base_forward <- probs_market(0,65-alder, increment_size, mu_base_forward)

for (i in 2:ncol(surrender_dat)){
  transprops_mu_base_forward$P00[,i] <- exp(-integrate_vec(surrender_dat[,i], mu_base = mu_base_forward))*ode_mu_base_forward[,"P00"]
  transprops_mu_base_forward$P01[,i] <- exp(-integrate_vec(surrender_dat[,i],  mu_base = mu_base_forward))*ode_mu_base_forward[,"P01"]
  transprops_mu_base_forward$P11[,i] <- exp(-integrate_vec(surrender_dat[,i],  mu_base = mu_base_forward))*ode_mu_base_forward[,"P11"]
  transprops_mu_base_forward$P10[,i] <- exp(-integrate_vec(surrender_dat[,i],  mu_base = mu_base_forward))*ode_mu_base_forward[,"P10"]
}


# mu_base = mu_as - true model --------------------------------------------

transprops_mu_base_true <- list(
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
    insert_time_first_col()
)

for (i in 2:ncol(surrender_dat)){
  mu_basetrue <-  approxfun(surrender_dat[,1], surrender_dat[,i]) 
  ode_mu_base_true <- probs_market(0,65-alder, increment_size, mu_basetrue)
  
  transprops_mu_base_true$P00[,i] <- ode_mu_base_true[,"P00"]
  transprops_mu_base_true$P01[,i] <- ode_mu_base_true[,"P01"]
  transprops_mu_base_true$P11[,i] <- ode_mu_base_true[,"P11"]
  transprops_mu_base_true$P10[,i] <- ode_mu_base_true[,"P10"]
}



# plots -------------------------------------------------------------------

ssh_plot_function <- 
  function(ssh,title , sti = NULL, diff = F, return_dat = F){
    
    time_grid <- seq(0, 65-alder, increment_size)
    
    if (is.null(sti)){
      
      true <- rowMeans(transprops_mu_base_true[[ssh]][,-1]) 
      forward <-rowMeans(transprops_mu_base_forward[[ssh]][,-1]) 
      mean <- rowMeans(transprops_mu_base_mean[[ssh]][,-1]) 
      simpel <- rowMeans(transprops_mu_base_0[[ssh]][,-1]) 
      
    }
    else{
      true <- transprops_mu_base_true[[ssh]][,sum(sti+1)]
      forward <-transprops_mu_base_forward[[ssh]][,sum(sti+1)]
      mean <- transprops_mu_base_mean[[ssh]][,sum(sti+1)]
      simpel <- transprops_mu_base_0[[ssh]][,sum(sti+1)]
    }
    
    if (diff){
      
      forward <- true-forward 
      mean <- true-mean
      simpel <- true-simpel
      
    }
    
    p <- ggplot(data.frame(time_grid),aes(time_grid)) +
      geom_line(aes(y=forward, colour="mu_base forward")) +
      geom_line(aes(y=mean, colour="mu_base mean")) +
      geom_line(aes(y=simpel, colour="mu_base 0")) +
      xlab("Tid") + ylab("Probability") + theme_bw()+
      ggtitle(title)
    
    if (!diff){
      p <- p+geom_line(aes(y=true, colour="mu_base true"))
    }
    
    if (return_dat)
      return(list(p=p, true=true, forward=forward, mean = mean ,simpel=simpel ))
    
    return(p)
  }

ssh_plot_function("P00", "P00 stig 1",6,diff = T)
grid.arrange(ssh_plot_function("P00", "",1,diff = T)  +  theme(legend.position = "none") + ylab("")+ xlab(""),
             ssh_plot_function("P01", "",1,diff = T) +  theme(legend.position = "none")+ ylab("")+ xlab(""),
             ssh_plot_function("P10", "",1,diff = T) +  theme(legend.position = "none")+ ylab("")+ xlab(""),
             ssh_plot_function("P11", "",1,diff = T) +  theme(legend.position = "none")+ ylab("")+ xlab(""),ncol=2,
             left = textGrob("Difference in path 1", rot = 90, vjust = 1),
             bottom = textGrob("Time", vjust = -1))

ssh_plot_function("P00", "P00", diff = F)
ssh_plot_function("P00", "P00 diff", diff = T)

ssh_plot_function("P01", "P01 stig 1", 2 ,diff = F)
ssh_plot_function("P01", "P01", diff = F)
ssh_plot_function("P01", "P01 diff", diff = T)

ssh_plot_function("P10", "P10 stig 1", 2 ,diff = F)
ssh_plot_function("P10", "P10", diff = F)
ssh_plot_function("P10", "P10 diff", diff = T)

ssh_plot_function("P11", "P11 stig 1", 2 ,diff = F)
ssh_plot_function("P11", "P11", diff = F)
ssh_plot_function("P11", "P11 diff", diff = T)

grid.arrange(ssh_plot_function("P00", "",diff = T)  +  theme(legend.position = "none")+ ylab("")+ xlab(""),
             ssh_plot_function("P01", "",diff = T) +  theme(legend.position = "none")+ ylab("")+ xlab(""),
             ssh_plot_function("P10", "",diff = T) +  theme(legend.position = "none")+ ylab("")+ xlab(""),
             ssh_plot_function("P11", "",diff = T) +  theme(legend.position = "none")+ ylab("")+ xlab(""),ncol=2,
             left = textGrob("Expected difference", rot = 90, vjust = 1),
             bottom = textGrob("Time", vjust = -1))


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

v_stjerne_bar_circ <- function(sti = NULL, diff=F){
  
  
  tek_reserve <- tek_reserve_cal(0, endtime-alder, increment_size) %>%
    filter(time<=65-alder)
  
  tmp_func <- function(ssh, sti){
    if (is.null(sti)){
      
      true <- rowMeans(transprops_mu_base_true[[ssh]][,-1]) 
      forward <-rowMeans(transprops_mu_base_forward[[ssh]][,-1]) 
      mean <- rowMeans(transprops_mu_base_mean[[ssh]][,-1]) 
      simpel <- rowMeans(transprops_mu_base_0[[ssh]][,-1]) 
      
    }
    else{
      true <- transprops_mu_base_true[[ssh]][,sum(sti+1)]
      forward <-transprops_mu_base_forward[[ssh]][,sum(sti+1)]
      mean <- transprops_mu_base_mean[[ssh]][,sum(sti+1)]
      simpel <- transprops_mu_base_0[[ssh]][,sum(sti+1)]
    }
    return(list(true=true, forward = forward, mean=mean, simpel=simpel))
  }
  
  P00 <- tmp_func("P00", sti)
  P01 <- tmp_func("P01", sti)
  P11 <- tmp_func("P11", sti)
  P10 <- tmp_func("P10", sti)
  
  V_z0_0_func <- function(char) P00[[char]]*rev(tek_reserve$V0)+P01[[char]]*rev(tek_reserve$V1)
  V_z0_1_func <- function(char) P11[[char]]*rev(tek_reserve$V1)+P10[[char]]*rev(tek_reserve$V0)
  
  #z0 =0
  V_z0_0_true <- V_z0_0_func("true")
  V_z0_0_forward <- V_z0_0_func("forward")
  V_z0_0_mean <- V_z0_0_func("mean")
  V_z0_0_simpel <- V_z0_0_func("simpel")
  
  #z0=1
  V_z0_1_true <- V_z0_1_func("true")
  V_z0_1_forward <- V_z0_1_func("forward")
  V_z0_1_mean <- V_z0_1_func("mean")
  V_z0_1_simpel <- V_z0_1_func("simpel")
  
  if(diff){
    out <-  cbind(time = rev(tek_reserve$time),V_z0_0_true- V_z0_0_forward,V_z0_0_true-V_z0_0_mean, 
                  V_z0_0_true-V_z0_0_simpel,
                  V_z0_1_true-V_z0_1_forward,V_z0_1_true-V_z0_1_mean,V_z0_1_true-V_z0_1_simpel ) %>%
      as.data.frame()
    colnames(out) <- c("time","V_z0_0_forward", "V_z0_0_mean", "V_z0_0_simpel", 
                       "V_z0_1_forward","V_z0_1_mean", "V_z0_1_simpel")
    return(out)
  }
  
  out <- cbind(time = rev(tek_reserve$time),V_z0_0_true, V_z0_0_forward,V_z0_0_mean, V_z0_0_simpel,
               V_z0_1_true,V_z0_1_forward,V_z0_1_mean,V_z0_1_simpel ) %>%
    as.data.frame()
  
  return(out)
  
}

v_stjerne_bar_circ_plot_function <- 
  function(sti = NULL, diff = F){
    
    time_grid <- seq(0, 65-alder, increment_size)
    
    tmp <- v_stjerne_bar_circ(sti, diff)
    
    p0 <- ggplot(data.frame(time_grid),aes(time_grid)) +
      geom_line(aes(y=tmp$V_z0_0_forward, colour="mu_base forward")) +
      geom_line(aes(y=tmp$V_z0_0_mean, colour="mu_base mean")) +
      geom_line(aes(y=tmp$V_z0_0_simpel, colour="mu_base 0")) +
      xlab("") + ylab("") + theme_bw()
    
    p1 <- ggplot(data.frame(time_grid),aes(time_grid)) +
      geom_line(aes(y=tmp$V_z0_1_forward, colour="mu_base forward")) +
      geom_line(aes(y=tmp$V_z0_1_mean, colour="mu_base mean")) +
      geom_line(aes(y=tmp$V_z0_1_simpel, colour="mu_base 0")) +
      xlab("") + ylab("") + theme_bw()
    
    if (!diff){
      p0 <- p0+geom_line(aes(y=tmp$V_z0_0_true, colour="mu_base true"))
      p1 <- p1+geom_line(aes(y=tmp$V_z0_1_true, colour="mu_base true"))
    }
    
    p <- grid.arrange(p0 +  theme(legend.position = "none"),
                      p1+ theme(legend.position = "none"),
                      left = textGrob("Expected difference", rot = 90, vjust = 1),
                      bottom = textGrob("Time", vjust = -1))
    
  }

v_stjerne_bar_circ_plot_function(diff=F)
v_stjerne_bar_circ_plot_function(diff=T)

v_stjerne_bar_circ_plot_function(1,diff=T)
v_stjerne_bar_circ_plot_function(1,diff=F)

