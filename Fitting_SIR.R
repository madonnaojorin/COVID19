##### Fitting the curve for different periods

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages ----------------------------------------------------------------
library(lubridate)
library(deSolve)
library(dplyr)

# Import case data
NAZ <- read.csv("Data/NAZ.csv")
CAZ <- read.csv("Data/CAZ.csv")
SAZ <- read.csv("Data/SAZ.csv")
NAZ$date <-  seq(as.Date(NAZ$date[1]), as.Date(NAZ$date[nrow(NAZ)]), "days")
CAZ$date <-  seq(as.Date(CAZ$date[1]), as.Date(CAZ$date[nrow(CAZ)]), "days")
SAZ$date <-  seq(as.Date(SAZ$date[1]), as.Date(SAZ$date[nrow(SAZ)]), "days")

# Import mobility data
NAZ_mob <- read.csv("Data/NAZ_mobility.csv")
CAZ_mob <- read.csv("Data/CAZ_mobility.csv")
SAZ_mob <- read.csv("Data/SAZ_mobility.csv")
NAZ_mob$date <-  seq(as.Date(NAZ_mob$date[1]), as.Date(NAZ_mob$date[nrow(NAZ_mob)]), "days")
CAZ_mob$date <-  seq(as.Date(CAZ_mob$date[1]), as.Date(CAZ_mob$date[nrow(CAZ_mob)]), "days")
SAZ_mob$date <-  seq(as.Date(SAZ_mob$date[1]), as.Date(SAZ_mob$date[nrow(SAZ_mob)]), "days")


save <- data.frame(sir_start_date=as.Date("2020-01-01"),sir_end_date=as.Date("2020-01-01"),
                   delay=10,WFH=0.5,
                   beta1=NA,beta2=NA,beta3=NA,
                   gamma1=NA,gamma2=NA,gamma3=NA,
                   convergence=NA,
                   value=NA)

DATE <- seq(as.Date("2020-03-01"),as.Date("2020-10-16"),"days")

# Set the initial parameter
(parameter0 = log(c(0.277347463,	0.864681623,	1.131945414,	0.197168164,	0.514870368,	0.751606605)))

for (i in 1:length(DATE)){
  # Set the period of interest 
  # Original: "2020-03-01" to "2020-10-17"
  sir_start_date <- DATE[i]
  sir_end_date <- "2020-10-17"
  WFH=0.5
  delay = 10
  
  mob_start_date <- as.Date(sir_start_date)-delay
  mob_end_date <- as.Date(sir_end_date)-delay
  
  mob1 <- subset(NAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$out_cases
  mob2 <- subset(CAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$out_cases
  mob3 <- subset(SAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$out_cases
  
  mob1_in <- subset(NAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$in_cases
  mob2_in <- subset(CAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$in_cases
  mob3_in <- subset(SAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$in_cases
  
  Mobility1 = as.data.frame(list(times = 1:(length(mob1)), mobility = mob1))
  input1 <- approxfun(Mobility1, rule = 2)
  Mobility2 = as.data.frame(list(times = 1:(length(mob2)), mobility = mob2))
  input2 <- approxfun(Mobility2, rule = 2)
  Mobility3 = as.data.frame(list(times = 1:(length(mob3)), mobility = mob3))
  input3 <- approxfun(Mobility3, rule = 2)
  
  Mobility1_in = as.data.frame(list(times = 1:(length(mob1_in)), mobility = mob1_in))
  input1_in <- approxfun(Mobility1_in, rule = 2)
  Mobility2_in = as.data.frame(list(times = 1:(length(mob2_in)), mobility = mob2_in))
  input2_in <- approxfun(Mobility2_in, rule = 2)
  Mobility3_in = as.data.frame(list(times = 1:(length(mob3_in)), mobility = mob3_in))
  input3_in <- approxfun(Mobility3_in, rule = 2)
  
  # 1: North Arizona 
  # 2: Central Arizona
  # 3: South Arizona
  pij <- read.csv("Data/commuteProp.csv")
  rownames(pij)<-pij$X
  p11 <- p22 <- p33 <- 1
  p12 <- pij["NAZ_R","CAZ_W"]*WFH
  p23 <- pij["CAZ_R","SAZ_W"]*WFH
  p32 <- pij["SAZ_R","CAZ_W"]*WFH
  p21 <- pij["CAZ_R","NAZ_W"]*WFH
  p13 <- pij["NAZ_R","SAZ_W"]*WFH
  p31 <- pij["SAZ_R","NAZ_W"]*WFH
  
  
  ##### SIR function 
  SIR_fun <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      parameters=exp(parameters)
      
      Mob1 <- 1/(1+exp(-input1(time)))
      Mob2 <- 1/(1+exp(-input2(time)))
      Mob3 <- 1/(1+exp(-input3(time)))
      
      Mob1_in <- 1/(1+exp(-input1_in(time)))
      Mob2_in <- 1/(1+exp(-input2_in(time)))
      Mob3_in <- 1/(1+exp(-input3_in(time)))
      
      foi1 <- (beta1*I1/N1*p11*Mob1+beta2*I2/N2*p12+beta3*I3/N3*p13)
      foi2 <- (beta1*I1/N1*p21+beta2*I2/N2*p22*Mob2+beta3*I3/N3*p23)
      foi3 <- (beta1*I1/N1*p31+beta2*I2/N2*p32+beta3*I3/N3*p33*Mob3)
      # Define transition equations
      # Population 1
      dS1 <- -S1*foi1
      dI1 <- S1*foi1-gamma1*I1
      dR1 <- gamma1*I1*Mob1_in
      
      # Population 2
      dS2 <- -S2*foi2
      dI2 <- S2*foi2-gamma2*I2
      dR2 <- gamma2*I2*Mob2_in
      
      # Population 3
      dS3 <- -S3*foi3
      dI3 <- S3*foi3-gamma3*I3
      dR3 <- gamma3*I3*Mob3_in
      list(c(dS1,dI1,dR1,dS2,dI2,dR2,dS3,dI3,dR3))
    })
  }
  
  
  # Put the daily cumulative incidence numbers into a vector called Infected
  Infected1 <- subset(NAZ, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$infected
  Infected2 <- subset(CAZ, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$infected
  Infected3 <- subset(SAZ, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$infected
  
  # Create an incrementing Day vector the same length as our cases vector
  Day <- 1:(length(Infected1))
  
  # Population of each region
  N1 = 212181 + 143476 + 	235099 + 21108 + 54018 + 110924 + 71887
  N2 = 4485414
  N3 = 213787 + 1047279 + 462789 + 38837 + 46498 + 9498 + 125922
  N  = N1+N2+N3
  
  # Specify initial values for S, I and R
  init <- c(
    S1 = N1 - Infected1[1],
    I1 = Infected1[1],
    R1 = 0,
    S2 = N2 - Infected2[1],
    I2 = Infected2[1],
    R2 = 0,  
    S3 = N3 - Infected3[1],
    I3 = Infected3[1],
    R3 = 0
  )
  
  # Define a function to calculate least squres
  MLE <- function(parameters) {
    
    names(parameters) <- c("beta1","beta2","beta3",
                           "gamma1","gamma2","gamma3")
    out <- ode(y = init, times = Day, func = SIR_fun, parms = exp(parameters),method = "bdf")
    out
    fit1 <- out[, "I1"]
    fit2 <- out[, "I2"]
    fit3 <- out[, "I3"]
    
    S1 <- sum((Infected1 - fit1)^2)
    S2 <- sum((Infected2 - fit2)^2)
    S3 <- sum((Infected3 - fit3)^2)
    S <- sum(S1,S2,S3)
    S
  }
  
  
  # Find the values of beta and gamma that give the
  # smallest MLE, which represents the best fit to the data.
  Opt <- try(optim(parameter0,
                   MLE,
                   method = "BFGS"))
  
  if(class(Opt)=="try-error"){
    Opt <- NA
    thisday <- max(save[save$convergence==0,]$sir_end_date)
    parameter0 <- log(save[save$sir_end_date==thisday,c("beta1","beta2","beta3","gamma1","gamma2","gamma3")])
  }else{
  # check for convergence
  Opt$convergence
  Opt_par <- setNames(exp(Opt$par), c("beta1","beta2","beta3",
                                      "gamma1","gamma2","gamma3"))
  Opt_par
  sir_start_date <- as.Date(sir_start_date)
  sir_end_date <- as.Date(sir_end_date)
  save[i,1] <- sir_start_date
  save[i,2] <- sir_end_date
  save[i,c(3,4)] <- c(delay,WFH)
  save[i,5:10] <- Opt_par[1:6]
  save[i,c(11,12)] <- c(Opt$convergence,Opt$value)
  parameter0 <- log(Opt_par[1:6])
  }
}

#write.csv(save,"results.csv")
