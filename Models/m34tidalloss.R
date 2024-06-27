m34tidalloss <- nimbleCode({
  # Priors (or fixed vals for things that may eventually be estimated)
  alpha ~ dunif(0,1)
  beta ~ dunif(1,10000)
  xinit ~ dunif(1,10000)
  theta ~ dunif(0.01,10)
  gamma ~ dnorm(0,10)
  
  
  # Process model
  
  # now add up expected DNA
  x[1] <- xinit
  for (t in 2:nTSo){
    # x is the expected (unobserved) concentration of DNA
    x[t] <- (x[t-1]+beta)*(alpha+gamma*tide_height_frac[t])
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
