m20lnormobservation <- nimbleCode({
  # Priors (or fixed vals for things that may eventually be estimated)
  alpha ~ dbeta(1, 5)
  beta ~ T(dlnorm(5.5, .5), 0, 10000) # truncated 
  xinit ~ T(dnorm(100, 50), 0, 1000) # truncated to not be negative
  theta ~ T(dnorm(0.025,1), 0, 10) # truncated to not be negative
  
  
  # Process model
  
  # now add up expected DNA
  x[1] <- xinit
  for (t in 2:nTSo){
    # x is the expected (unobserved) concentration of DNA
    x[t] <- (x[t-1]+beta)*alpha
  }
  
  # Observation model
  for (t in 1:nTSo){
    # calculate sd given CV
    sigma[t] <- theta*x[t]
    sd_log_x[t] <- sqrt(log(1+(sigma[t]^2/x[t]^2)))
    # y are the observed concentrations of DNA
    # y[t] ~ dnorm(mean=x[t],sd=sigma[t])
    y[t] ~ dlnorm(mean=log(x[t]),sd=sd_log_x[t])
  }
  
})
