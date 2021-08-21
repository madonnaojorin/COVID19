SEIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    Mob1 <- 1/(1+exp(-input1(time)))
    Mob2 <- 1/(1+exp(-input2(time)))
    Mob3 <- 1/(1+exp(-input3(time)))
    #print(c(Mob1,Mob2,Mob3))
    
    nu1 <- NAZ_vaccine_rate(time) 
    nu2 <- CAZ_vaccine_rate(time)
    nu3 <- SAZ_vaccine_rate(time)
    # print(c(nu1,nu2,nu3))
    
    foi1 <- (beta1*I1/N1*p11*Mob1+beta2*I2/N2*p12+beta3*I3/N3*p13)
    foi2 <- (beta1*I1/N1*p21+beta2*I2/N2*p22*Mob2+beta3*I3/N3*p23)
    foi3 <- (beta1*I1/N1*p31+beta2*I2/N2*p32+beta3*I3/N3*p33*Mob3)
    
    # Define transition equations
    # Population 1
    dS1 <- -S1 * foi1 - nu1 * S1
    dE1 <- S1 * foi1 - sigma1 * E1
    dI1 <- sigma1 * E1 - gamma1 * I1
    dR1 <- gamma1 * I1 + nu1 * S1
    
    # Population 2
    dS2 <- -S2 * foi2 - nu2 * S2
    dE2 <- S2 * foi2 - sigma2 * E2
    dI2 <- sigma2 * E2 - gamma2 * I2
    dR2 <- gamma2 * I2 + nu2 * S2
    
    # Population 3
    dS3 <- -S3 * foi3 - nu3 * S3
    dE3 <- S3 * foi3 - sigma3 * E3
    dI3 <- sigma3 * E3 - gamma3 * I3
    dR3 <- gamma3 * I3 + nu3 * S3
    
    list(c(dS1,dE1,dI1,dR1,dS2,dE2,dI2,dR2,dS3,dE3,dI3,dR3))
  })
}