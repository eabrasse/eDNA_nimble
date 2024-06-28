m35randomloss <- nimbleCode({
  # Priors (or fixed vals for things that may eventually be estimated)
  beta ~ dunif(1,10000)
  xinit ~ dunif(1,10000)
  theta ~ dunif(0.01,10)
  
  
  # Process model
  
  # now add up expected DNA
  x[1] <- xinit
  for (t in 2:nTSo){
    # x is the expected (unobserved) concentration of DNA
    x[t] <- alpha[t]*x[t-1]+beta
  }
  
  # Observation model
  for (t in 1:nTSo){
    # calculate sd given CV
    sigma[t] <- theta*x[t]
    sd_log_x[t] <- sqrt(log(1+(sigma[t]^2/x[t]^2)))
    y[t] ~ dlnorm(mean=log(x[t]),sd=sd_log_x[t])
    alpha[t] ~ dunif(0,1)
  }
  
})
