# by Hannah Meredith
# Updated: October 24, 2018

# libraries

# install.packages("pracma")
library("pracma")
# install.packages("deSolve")
library(deSolve)
# install.packages("ggplot2")
library("ggplot2")
# install.packages("readxl")
library(readxl)
library(grid)
library(gridExtra)
library(reshape2)


# baseline (no controls)
baseline <- function(y0, P, m, p_h, p_n) {
  P = c(
    P,
    C_n = 0 ,
    C_l = 0 ,
    C_h = 0 ,
    m = m,
    p_h = p_h,
    dose_delay = 0
  )
  
  step = 1
  
  t = seq(0, P[13] + P[14], by = step) #duration of simulation = t_ss + treatment period
  
  base <-
    rk(
      y = y0,
      times = t,
      func = net_shift_biting,
      parms = P,
      method = "ode45",
      atol = 1e-10,
      rtol = 1e-10
    )
  
  total_infected_h <- base[, 3] + base[, 5]
  time <- base [, 1]
  baseline <- list(time, total_infected_h)
  
}

# bednets alone
nets <- function(y0, P, m, p_h, p_n) {
  P = c(
    P,
    C_n = 1 * 0.75 ,
    C_l = 0 ,
    C_h = 0 ,
    period_D = 365,
    m = m,
    p_h = p_h,
    dose_delay = 0
  )
  
  nets <- dosing(y0, P)
}

# Dosing function

dosing <- function(y0, P) {
  
  t_steadyState <- as.vector(P[13])
  treatment_pd <- as.vector(P[14])
  period_N <- P[16]
  net <- P[17]
  dose_l <- P[26]
  dose_h <- P[27]
  period_D <- P[33]
  dose_delay <- P [36]
  
  # Initialize time window to achieve steady state and run ODEs
  step <- 1
  times = seq(0, t_steadyState, by = step)
  out <-
    rk(
      y = y0,
      times = times,
      func = net_shift_biting,
      parms = P,
      method = "ode45",
      atol = 1e-10,
      rtol = 1e-10
    )
  overall <- out
  
  # After system reaches steady state, introduce control methods
  
  # Define new time windows. These start and end times for nets and drug dosing
  # will determine the order of subsequent doses
  t_startD = as.vector(tail(overall[, 1], n = 1) + dose_delay)
  t_endD = as.vector(tail(overall[, 1], n = 1) + dose_delay + period_D)
  t_startN = as.vector(tail(overall[, 1], n = 1))
  t_endN = as.vector(tail(overall[, 1], n = 1) + period_N)
  
  # Update counters
  doseNum = 0
  netNum = 0
  toggle = 0
  
  # Update initial conditions for next ODE solver. If no drug delay, add first nets and drug; otherwise just add nety)
  if (dose_delay == 0){
    y0 <- c(
      E_h = as.vector(tail(overall[, 2], n = 1)),
      I_h = as.vector(tail(overall[, 3], n = 1)),
      R_h = as.vector(tail(overall[, 4], n = 1)),
      A_h = as.vector(tail(overall[, 5], n = 1)),
      E_m = as.vector(tail(overall[, 6], n = 1)),
      I_m = as.vector(tail(overall[, 7], n = 1)),
      D_l = as.vector(tail(overall[, 8], n = 1) + dose_l),
      D_h = as.vector(tail(overall[, 9], n = 1) + dose_h),
      N = as.vector(net))
    
    out <-
      rk(
        y = y0,
        times = seq(t_startD, t_endD, by = step),
        func = net_shift_biting,
        parms = P,
        method = "ode45",
        atol = 1e-10,
        rtol = 1e-10
      )
    doseNum = doseNum + 1
    netNum = netNum + 1
  } else {
    y0 <- c(
      E_h = as.vector(tail(overall[, 2], n = 1)),
      I_h = as.vector(tail(overall[, 3], n = 1)),
      R_h = as.vector(tail(overall[, 4], n = 1)),
      A_h = as.vector(tail(overall[, 5], n = 1)),
      E_m = as.vector(tail(overall[, 6], n = 1)),
      I_m = as.vector(tail(overall[, 7], n = 1)),
      D_l = as.vector(tail(overall[, 8], n = 1)),
      D_h = as.vector(tail(overall[, 9], n = 1)),
      N = as.vector(net))
    
    out <-
      rk(
        y = y0,
        times = seq(t_startN, t_startD, by = step),
        func = net_shift_biting,
        parms = P,
        method = "ode45",
        atol = 1e-10,
        rtol = 1e-10
      )
    netNum = netNum + 1
  }
  
  #Add new times to overall timecourse matrix
  overall <- rbind(overall[1:nrow(overall)-1,], out)
  t_startD = as.vector(tail(overall[, 1], n = 1))
  t_endD = as.vector(tail(overall[, 1], n = 1) + period_D)
  
  while (tail(overall[, 1], n = 1) < t_steadyState + treatment_pd &
         toggle == 0) {
    # Ensures simulation cuts off at 10 years
    if (t_endD >= t_steadyState + treatment_pd) {
      t_endD = t_steadyState + treatment_pd
      t_endN = t_steadyState + treatment_pd
      toggle = 1
    }
    y0 <- c(
      E_h = as.vector(tail(overall[, 2], n = 1)),
      I_h = as.vector(tail(overall[, 3], n = 1)),
      R_h = as.vector(tail(overall[, 4], n = 1)),
      A_h = as.vector(tail(overall[, 5], n = 1)),
      E_m = as.vector(tail(overall[, 6], n = 1)),
      I_m = as.vector(tail(overall[, 7], n = 1)),
      D_l = as.vector(tail(overall[, 8], n = 1) + dose_l),
      D_h = as.vector(tail(overall[, 9], n = 1) + dose_h),
      N = as.vector(tail(overall[, 10], n = 1))
    )
    out <-
      rk(
        y = y0,
        times = seq(t_startD, t_endD, by = step),
        func = net_shift_biting,
        parms = P,
        method = "ode45",
        atol = 1e-10,
        rtol = 1e-10
      )
    overall <- rbind(overall[1:nrow(overall)-1,], out)
    t_startD = as.vector(tail(overall[, 1], n = 1))
    t_endD = as.vector(tail(overall[, 1], n = 1) + period_D)
    doseNum = doseNum + 1
  }
  
  TI_h = overall[, 3] + overall[, 5]
  overall_TI_h <- cbind(overall, TI_h)
  overall_TI_h.df <- as.data.frame(overall_TI_h)
  #
  ggplot(overall_TI_h.df, aes(x = time)) +
    geom_line(aes(y = TI_h, colour = "Total infected")) +
    geom_line(aes(y = E_h, colour = "Exposed")) +
    geom_line(aes(y = I_h, colour = "Infected")) +
    geom_line(aes(y = R_h, colour = "Recovered")) +
    geom_line(aes(y = A_h, colour = "Asymptomatically Infected")) +
    ylab(label = "% Human hosts") +
    xlab(label = "Time") +
    coord_cartesian(xlim = c(1000, 1500), ylim = c(0, 1))
  #
  ggplot(overall_TI_h.df, aes(x = time)) +
    geom_line(aes(y = D_l, colour = "Drug")) +
    # geom_line(aes(y = N, colour = "Net")) +
    ylab(label = "Drug") +
    xlab(label = "Time") +
    coord_cartesian(xlim = c(0, 5000), ylim = c(0, 10))
  
  total_infected_h <- overall[, 3] + overall[, 5]
  time <- overall [, 1]
  dosing <- list(time, total_infected_h, overall[, 10])
}

# ODEs for malaria transmission with control methods (LLINs and Systemic insecticides)
net_shift_biting <- function(t, y, parms) {
  with(as.list(c(parms, y)), {
    k_n = log(2) / hl_n
    k_d_l = log(2) / hl_d_l
    k_d_h = log(2) / hl_d_h
    
    
    # EQs and ODEs
    phi <- N_t/(N_t + N)
    b_h <-
      a * b * phi * p_h * (1 -  phi * p_n * C_n * N ^ hill_n / (N ^ hill_n + N_death ^ hill_n))
    b_m <-
      a * c * phi * p_h * (1 - phi * p_n * C_n * N ^ hill_n / (N ^ hill_n + N_death ^ hill_n))
    mu_2c <-
      mu_2 + 1 / 3 * (
        phi * p_h * phi * p_n * C_n * N ^ hill_n / (N ^ hill_n + N_death ^ hill_n) * mu_n + 
          (1 - phi* p_h) * C_l * D_l ^hill_d / (D_l ^ hill_d + D_death ^ hill_d) *mu_d + 
          (phi * p_h) * C_h * D_h ^ hill_d / (D_h ^ hill_d + D_death ^hill_d) * mu_d)
    
    
    S_h <- 1 - E_h - I_h - R_h - A_h
    S_m <- 1 - E_m - I_m
    
    dE_h <-
      m * b_h * S_h * I_m  - (1 / tau_h + r + mu_1) * E_h
    dI_h <- 1 / tau_h * E_h - (r + q1 + mu_1) * I_h
    dR_h <-
      q1 * (I_h + A_h) - (theta * b_h * m  * I_m + q2 + mu_1) * R_h
    dA_h <-
      theta * b_h * m * I_m * R_h - (q1 + r + mu_1) * A_h
    
    dE_m <-
      b_m * (I_h + sigma * A_h) * S_m - (1 / tau_m + mu_2c) * E_m
    dI_m <- (1 / tau_m) * E_m - mu_2c * I_m
    dD_l <- -k_d_l * D_l
    dD_h <- -k_d_h * D_h
    dN <- -k_n * N
    
    # print(c(t,mu_2c))
    
    res <- c(dE_h, dI_h, dR_h, dA_h, dE_m, dI_m, dD_l, dD_h, dN)
    list(res)
  })
}

# define parameters and  ranges

# Parameters
P <- c(
  a = 0.2,  #1
  b = 0.5,  #2
  c = 0.5,  #3
  r = 0.01, #4
  mu_1 = 1 / 21900,  #5
  mu_2 = 0.12,  #6
  tau_m = 10,   #7
  tau_h = 21,   #8
  q1 = 1 / 200, #9
  q2 = 1 / 1000,#10 
  sigma = 0.25, #11
  theta = 0.5,  #12
  t_ss = 1000, # 13 time to let system reach steady state
  treatment_pd = 4 * 365, # 14 period over which the treatment is analyzed
  N_death = 0.73, # 15 LC50 of LLIN for resistant mosquitoes
  period_N = 1095, # 16 net is replaced every 3 years
  net = 2, # 17 dose of new net
  hl_n = 1906, # 18 halflife of bednet insecticide
  mu_d = 0.6, # 19 death rate due to drug
  mu_n = 0.1, # 20 death rate due to bed net
  hill_n = 2, # 21
  hill_d = 2, # 22
  N_t = 1.34, # 23 threshold LLIN concentration at which behaviour returns to pre-net status
  hl_d_l = 5.1,  # 24 half life of ivermectin subcutaneously applied to cattle
  hl_d_h = 3.4,  # 25 half life of ivermectin orally applied to humans
  dose_l = 28.1, # 26 cmax of ivermectin subcutaneously applied to cattle
  dose_h = 43.19, # 27 cmax of ivermectin subcutaneously applied to human
  D_death = 7.5,  # 28 LC50 of Anopheles to ivermectin,
  p_n = 0.2 #29 
)

period_D_range <- c(365/12, 365) # dosing frequencies to compare
period_title <- c("Monthly dosing", "Yearly dosing")
dose_delay_range <- c(0, 365/12, 365, 2*365, 3*365) #dose delays to compare

# define initial conditions
y0 <- c(
  E_h = 0.01,
  I_h = 0.01,
  R_h = 0.01,
  A_h = 0.01,
  E_m = 0.01,
  I_m = 0.01,
  D_l = 0,
  D_h = 0,
  N = 0
)

# calculate number of infected cases without any control methods (baseline)
base <- baseline(y0, P, m = 20, p_h = 0.8, p_n = 0.2)
time_base <- base[[1]]
infected_base <- base[[2]]
infected_base_steadyState <- tail(base[[2]], n = 1)

# calculate number of infected cases with only bed nets
net <- nets(y0, P, m = 20, p_h = 0.8, p_n = 0.2)
time_net <- net[[1]]
infected_net <- net[[2]]
net_dosing <- net[[3]]
net.df <- data.frame(time_net, infected_net)
infected_net_steadyState <- tail(net[[2]], n = 1)
prevalence_ratio_net <-
  infected_net / infected_base_steadyState # used in plotting time courses

pltlist <- list()

# calculate prevalence with bednets and different systemic insecticide treatments
for (i in size(period_D_range, 2):1) {
  x1 <- time_net
  y1 <- prevalence_ratio_net
  df1 <- data.frame(x1, y1)
  df1["Delay"] <- "LLIN only"
  
  for (j in size(dose_delay_range,2):1) {
  P = c(
    P[1:29],
    C_n = 1 * 0.75,
    C_l = 1 * 1,
    C_h = 1 * 1,
    period_D = period_D_range[i],
    m = 20,
    p_h = 0.8,
    dose_delay = dose_delay_range[j]
   )
  
  # calculate number of infected cases with bednets and systemic insecticides
  # compute the relative prevalence of malaria cases of avoided with addition of SI
  all <- dosing(y0, P)
  time_all <- all[[1]]
  infected_all <- all[[2]]
  net_dosing <- all[[3]]
  
  all.df <- data.frame(time_all, infected_all)
  
  prevalence_ratio_all <-
    infected_all / infected_base_steadyState # used in plotting time courses
  
  x1 <- time_all
  y1 <- prevalence_ratio_all
  
  df2 <- data.frame(x1,y1)
  df2["Delay"] <- toString(round(dose_delay_range[j]))
  df1 <-rbind(df1,df2)
 
  
  df_year<-df1
  df_year$x1 <- (df1$x1 - P[13])/365
  }
  
  p <- ggplot(data = df_year) + geom_line(aes(x1, y1, colour = Delay), size=1) + 
    theme(aspect.ratio = 1) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    xlab("Time (years)") +
    theme(
      axis.text.x = element_text(colour = 'black', size = 10),
      axis.title.x = element_text(
        size = 12,
        hjust = 0.5,
        vjust = 0.5
      )
    ) +
    ylab("Malaria prevalence ratio") +
    theme(
      axis.text.y = element_text(colour = 'black', size = 10),
      axis.title.y = element_text(
        size = 12,
        hjust = 0.5,
        vjust = 0.2
      )
    ) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 1), expand=FALSE)+
    ggtitle(toString((period_title[i])))+
    theme(plot.title = element_text(size = 12))
  
  pltlist[[i]] <- (p) 
}

# get legend for common plot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(pltlist[[1]])


# subplots
grid.arrange(
  pltlist[[1]] + theme(legend.position = "none"),
  pltlist[[2]] + theme(legend.position = "none"),
  legend,
  nrow = 1,
  ncol = 3
)

