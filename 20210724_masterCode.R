##### Master code:
# Created on July 24 2021, modified on August 5 2021
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages ----------------------------------------------------------------
libraries = c("magrittr","dplyr","ggplot2","lubridate","deSolve","gridExtra","scales")
for (x in libraries) {
  library(x, character.only = TRUE, warn.conflicts = FALSE)
}

# Plot settings for all the figures----------------------------------------
Font1 = 40; Font2 = 38; Font3 = 28

###### Set the date to analyze
sir_start_date <- "2020-09-12"
sir_end_date <- "2021-06-30"
delay = 10
mob_start_date <- as.Date(sir_start_date)-delay
mob_end_date <- as.Date(sir_end_date)-delay

##### Import data 
# COVID-19 case data
NAZ <- read.csv("Data/Final_data/NAZ.csv")
CAZ <- read.csv("Data/Final_data/CAZ.csv")
SAZ <- read.csv("Data/Final_data/SAZ.csv")
NAZ$date <-  seq(as.Date(NAZ$date[1]), as.Date(NAZ$date[nrow(NAZ)]), "days")
CAZ$date <-  seq(as.Date(CAZ$date[1]), as.Date(CAZ$date[nrow(CAZ)]), "days")
SAZ$date <-  seq(as.Date(SAZ$date[1]), as.Date(SAZ$date[nrow(SAZ)]), "days")
print(c(NAZ$date[1],CAZ$date[1],SAZ$date[1],NAZ$date[nrow(NAZ)],CAZ$date[nrow(CAZ)],SAZ$date[nrow(SAZ)]))

# Subset the data of targeted period
Infected1 <- subset(NAZ, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$infected
Infected2 <- subset(CAZ, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$infected
Infected3 <- subset(SAZ, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$infected

######  Extract the mobility data ###### 
NAZ_mob <- read.csv("Data/Final_data/NAZ_mobility.csv")
CAZ_mob <- read.csv("Data/Final_data/CAZ_mobility.csv")
SAZ_mob <- read.csv("Data/Final_data/SAZ_mobility.csv")
NAZ_mob$date <-  seq(as.Date(NAZ_mob$date[1]), as.Date(NAZ_mob$date[nrow(NAZ_mob)]), "days")
CAZ_mob$date <-  seq(as.Date(CAZ_mob$date[1]), as.Date(CAZ_mob$date[nrow(CAZ_mob)]), "days")
SAZ_mob$date <-  seq(as.Date(SAZ_mob$date[1]), as.Date(SAZ_mob$date[nrow(SAZ_mob)]), "days")
print(c(NAZ_mob$date[1],CAZ_mob$date[1],SAZ_mob$date[1],NAZ_mob$date[nrow(NAZ_mob)],CAZ_mob$date[nrow(CAZ_mob)],SAZ_mob$date[nrow(SAZ_mob)]))

# Data stored as negative/positive. Convert it to percentage 
(mob1 <- 1 + subset(NAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$out_cases)
(mob2 <- 1 + subset(CAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$out_cases)
(mob3 <- 1 + subset(SAZ_mob, date >= ymd(mob_start_date) & date <= ymd(mob_end_date))$out_cases)

# Interpolate mobility data
Mobility1 = as.data.frame(list(times = 1:(length(mob1)), mobility = mob1))
input1 <- approxfun(Mobility1, rule = 2)
Mobility2 = as.data.frame(list(times = 1:(length(mob2)), mobility = mob2))
input2 <- approxfun(Mobility2, rule = 2)
Mobility3 = as.data.frame(list(times = 1:(length(mob3)), mobility = mob3))
input3 <- approxfun(Mobility3, rule = 2)

###### Extract the vaccine data ######
NAZ_vaccine_df <- read.csv("Data/Final_data/NAZ_vaccine.csv")[,-1]
CAZ_vaccine_df <- read.csv("Data/Final_data/CAZ_vaccine.csv")[,-1]
SAZ_vaccine_df <- read.csv("Data/Final_data/SAZ_vaccine.csv")[,-1]

vaccine20 <- data.frame(date = seq(as.Date(NAZ$date[1]), as.Date("2020-12-12"), "days"),
                        population = 0,
                        first_dose = 0,
                        first_dose_prop = 0,
                        first_dose_cum = 0,
                        fully = 0,
                        fully_prop = 0,	
                        fully_cum = 0,
                        group = rep("Vac20")
)

NAZ_vaccine <- rbind(vaccine20,NAZ_vaccine_df)%>%mutate(zero = 0)
CAZ_vaccine <- rbind(vaccine20,CAZ_vaccine_df)%>%mutate(zero = 0)
SAZ_vaccine <- rbind(vaccine20,SAZ_vaccine_df)%>%mutate(zero = 0)
print(c(NAZ_vaccine$date[1],CAZ_vaccine$date[1],SAZ_vaccine$date[1],NAZ_vaccine$date[nrow(NAZ_vaccine)],CAZ_vaccine$date[nrow(CAZ_vaccine)],SAZ_vaccine$date[nrow(SAZ_vaccine)]))

# Interporate mobility data
# Use fully_prop which is the vaccine rate calculated by fully vaccinated population/population/day
(NAZ_vac <- subset(NAZ_vaccine, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$fully_prop)
NAZ_vaccine_rate <- approxfun(NAZ_vac, rule = 2)
(CAZ_vac <- subset(CAZ_vaccine, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$fully_prop)
CAZ_vaccine_rate <- approxfun(CAZ_vac, rule = 2)
(SAZ_vac <- subset(SAZ_vaccine, date >= ymd(sir_start_date) & date <= ymd(sir_end_date))$fully_prop)
SAZ_vaccine_rate <- approxfun(SAZ_vac, rule = 2)

###### Extract the Commute data
# 1: North Arizona, 2: Central Arizona, 3: South Arizona
pij <- read.csv("Data/Final_data/commuteProp.csv")
# Multiplied by WFH factor calculated from residential data (WFH_calculation.R)
rownames(pij)<-pij$X
#WFH=0.5
p11 <- p22 <- p33 <- 1
p12 <- pij["NAZ_R","CAZ_W"]*(1-0.05739594)
p23 <- pij["CAZ_R","SAZ_W"]*(1-0.0881295)
p32 <- pij["SAZ_R","CAZ_W"]*(1-0.06303537)
p21 <- pij["CAZ_R","NAZ_W"]*(1-0.0881295)
p13 <- pij["NAZ_R","SAZ_W"]*(1-0.05739594)
p31 <- pij["SAZ_R","NAZ_W"]*(1-0.06303537)

# Create a Day vector the same length as our cases vector
Day <- 1:(length(Infected1))

# Calculate the total population of each region 
N1 = 219583 + 145326 + 236426 + 22542 + 55293 + 113276 + 72644 
N2 = 4485414
N3 = 235321 + 1052375 + 467932 + 38666 + 53731 + 10558 + 131694 
N  = N1+N2+N3

# Specify initial values for N, S, E, I and R
# Susceptible population is total population minus cumulative number of cases on 2020-09-11
(init <- c(
  S1 = N1 - Infected1[1] - NAZ$cum_cases[NAZ$date == "2020-09-11"],
  E1 = 0,
  I1 = Infected1[1],
  R1 = NAZ$cum_cases[NAZ$date == "2020-09-11"],
  
  S2 = N2 - Infected2[1] - CAZ$cum_cases[CAZ$date == "2020-09-11"],
  E2 = 0,
  I2 = Infected2[1],
  R2 = CAZ$cum_cases[CAZ$date == "2020-09-11"],  
  
  S3 = N3 - Infected3[1] - SAZ$cum_cases[SAZ$date == "2020-09-11"],
  E3 = 0,
  I3 = Infected3[1],
  R3 = SAZ$cum_cases[SAZ$date == "2020-09-11"]
))

# Set the initial values for parameters 
(parameter0= c(0.803, 0.116, 0.376, 0.906, 0.110, 0.377, 0.289, 0.925, 0.216)) 
parameter0 <- setNames(parameter0, c("beta1","sigma1","gamma1",
                                     "beta2","sigma2","gamma2",
                                     "beta3","sigma3","gamma3"))

# Load function files for SEIR and MLE
source(file = "Script/RSS_function.R")
source(file = "Script/SEIR_function.R")

Opt <- optim(parameter0,
             RSS,
             method = "L-BFGS-B",lower = rep(0, length(parameter0)), upper = rep(1, length(parameter0)),control = list(maxit = 10000))

# check for convergence
Opt$convergence

Opt_par <- setNames(Opt$par, c("beta1","sigma1","gamma1",
                               "beta2","sigma2","gamma2",
                               "beta3","sigma3","gamma3"))
Opt_par

save <- c(Opt_par[1],Opt_par[2],Opt_par[3],
          Opt_par[4],Opt_par[5],Opt_par[6],
          Opt_par[7],Opt_par[8],Opt_par[9],
          Opt$convergence,Opt$value)

write.csv(save,"Results/Results0727_mac.csv")

(parameters <- Opt_par[1:9])
t <- Day

fitted_prevalence <- data.frame(ode(
  y = init, times = t,
  func = SEIR, parms = parameters
))

################## NAZ ##################
fitted_prevalence <- fitted_prevalence %>%
  mutate(
    Date = ymd(sir_start_date) + days(t - 1),
    Country = "NAZ",
    Prevalence = Infected1
  )

data2plot <- data.frame(Date = as.Date(rep(fitted_prevalence$Date,2)),
                        data = c(fitted_prevalence$I1,fitted_prevalence$Prevalence),
                        group = c(rep("Estimated",nrow(fitted_prevalence)),rep("Observed",nrow(fitted_prevalence))))

NAZ_pt <- ggplot(data2plot, aes(x = Date, y=data, group = group)) +
  geom_line(aes(linetype=group, color = group), size = 1.3) +
  labs(y = "Prevalence",
       title ="NAZ")+
  scale_y_continuous(limits=c(0,100000),labels = comma)+
  scale_x_date(labels = NULL, breaks = NULL)+
  labs(x = "")+
  # Delete x for figures together
  # scale_x_date(breaks = "2 week",
  #              labels = date_format("%m/%d"),
  #              limits = c(as.Date(sir_start_date),as.Date(sir_end_date))) +
  scale_color_manual(values=c("#00468BFF","#ED0000FF"))+
  theme(
    panel.background = element_blank(),
    plot.title = element_text(face = "bold",size=Font1),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size =Font3),
    axis.title = element_text(size =Font2, face = "bold"),
    axis.text.x = element_text(size = Font3, angle = 45, vjust = 0.5), # Delete for 3 figures
    legend.key = element_blank(),
    legend.text = element_text(size=Font3),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.position = "NONE",
    legend.title = element_blank()) 
NAZ_pt
################## CAZ ##################
fitted_prevalence2 <- fitted_prevalence %>%
  mutate(
    Date = ymd(sir_start_date) + days(t - 1),
    Country = "CAZ",
    Prevalence = Infected2
  )

# plot the data
data2plot2 <- data.frame(Date = rep(fitted_prevalence2$Date,2),
                         data = c(fitted_prevalence2$I2,fitted_prevalence2$Prevalence),
                         group = c(rep("Estimated",nrow(fitted_prevalence2)),rep("Observed",nrow(fitted_prevalence2))))


CAZ_pt <- ggplot(data2plot2, aes(x = Date, y=data, group = group)) +
  geom_line(aes(linetype=group, color = group), size = 1.3) +
  labs(y = "Prevalence",
       title = "CAZ")+
  scale_y_continuous(limits=c(0,100000),labels = comma)+
  scale_x_date(labels = NULL, breaks = NULL)+
  labs(x = "")+
  # Delete x for figures together
  # scale_x_date(breaks = "2 week",
  #              labels = date_format("%m/%d"),
  #              limits = c(as.Date(sir_start_date),as.Date(sir_end_date))) +
  scale_color_manual(values=c("#00468BFF","#ED0000FF"))+
  theme(
    panel.background = element_blank(),
    plot.title = element_text(face = "bold",size=Font1),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size =Font3),
    axis.title = element_text(size =Font2, face = "bold"),
    axis.text.x = element_text(size = Font3, angle = 45, vjust = 0.5), # Delete for 3 figures
    legend.key = element_blank(),
    legend.text = element_text(size=Font3),
    legend.spacing.x = unit(1.0, 'cm'),
    legend.position = "NONE",
    legend.title = element_blank()) 
CAZ_pt

################## SAZ ##################
fitted_prevalence3 <- fitted_prevalence %>%
  mutate(
    Date = ymd(sir_start_date) + days(t - 1),
    Country = "SAZ",
    Prevalence = Infected3
  )

# plot the data
data2plot3 <- data.frame(Date = rep(fitted_prevalence3$Date,2),
                         data = c(fitted_prevalence3$I3,fitted_prevalence3$Prevalence),
                         group = c(rep("Estimated",nrow(fitted_prevalence3)),rep("Observed",nrow(fitted_prevalence3))))


SAZ_pt <- ggplot(data2plot3, aes(x = Date, y=data, group = group)) +
  geom_line(aes(linetype=group, color = group), size = 1.3) +
  labs(y = "Prevalence",
       title = "SAZ")+
  scale_y_continuous(limits=c(0,100000),labels = comma)+
  scale_color_manual(values=c("#00468BFF","#ED0000FF"))+
  scale_x_date(breaks = "2 week",
               labels = date_format("%m/%d"),
               limits = c(as.Date(sir_start_date),as.Date(sir_end_date))) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size =Font3),
    axis.title = element_text(size =Font2, face = "bold"),
    axis.text.x = element_text(size = Font3, angle = 45, vjust = 0.5),
    plot.title = element_text(face = "bold",size=Font1),
    legend.key = element_blank(),
    legend.text = element_text(size=Font3),
    legend.position = "NONE",
    legend.title = element_blank()) 
SAZ_pt

Three <- grid.arrange(NAZ_pt,CAZ_pt,SAZ_pt,ncol = 1)
ggsave(filename = "Results/Result0726.png", Three, dpi = 400,width = 300,height = 600,units = c("mm"))


