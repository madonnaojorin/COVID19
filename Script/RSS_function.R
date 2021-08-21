RSS <- function(parameters) {
  
  names(parameters) <- c("beta1","sigma1","gamma1",
                         "beta2","sigma2","gamma2",
                         "beta3","sigma3","gamma3")
  out <- ode(y = init, times = Day, func = SEIR, parms = parameters, method = "ode45") 
  
  fit1 <- out[, "I1"]
  fit2 <- out[, "I2"]
  fit3 <- out[, "I3"]
  
  S1 <- sum((Infected1 - fit1)^2)
  S2 <- sum((Infected2 - fit2)^2)
  S3 <- sum((Infected3 - fit3)^2)
  #print(paste("S1:", S1))
  #print(paste("S2:", S2))
  #print(paste("S3:", S3))
  S <- sum(S1,S2,S3)
  S
  print(S)
}