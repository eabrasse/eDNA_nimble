

library(nimble)

source("./Models/m21lineardolphin.R")

moddata <- read.csv("../../../data/dolphin_tide_minutefreq.csv")
# moddata <- read.csv("../../../data/ESP_box_model_terms_ESPtimesteps.csv")
obsdata <- read.csv("../../../data/ESP_timestamps.csv")

nimbleData <- list(
  y = obsdata$PB_quantity_mean
  
)

t0 <- c()
t1 <- c()
N_dolphins <- c()
for (t in 1:nrow(obsdata)) {
  t0[t] <- min(which(moddata$timestamp > obsdata$timestamp_0[t]))
  t1[t] <- max(which(moddata$timestamp < obsdata$timestamp_1[t]))
  N_dolphins[t] <- mean(moddata$ndolphin_net[t0[t]:t1[t]])
}


nimbleConsts <- list(
  #t0 = t0,
  #t1 = t1,
  #nTSm = nrow(moddata),
  nTSo = nrow(obsdata),
  N_dolphins = N_dolphins
  #Mod_ts = moddata$timestamp
  #theta = 0.27
)

nimbleInits <- list(alpha = 0.2,
                    beta = 75,
                    xinit = 100,
                    x = rep(1, nrow(obsdata)),
                    sigma = rep(1, nrow(obsdata)),
                    theta = 1.0
)

nimbleParams <- list("x", "alpha", "beta", "xinit", "sigma","sd_log_x","theta")

model <- nimbleModel(code = m21lineardolphin,
                     data = nimbleData,
                     constants = nimbleConsts,
                     inits = nimbleInits,
                     check = FALSE)

set.seed(20240401) 

nimbleOut <- nimbleMCMC(model, 
                        monitors = nimbleParams, 
                        thin = 100, niter = 1000000, nburnin = 750000, nchains = 4,
                        summary = TRUE,WAIC = TRUE)

save(nimbleOut, file = "./Results/nimbleOut_m21lineardolphin.RData")

# load("./Results/nimbleOut_m21lineardolphin.RData")

nimbleAC <- as.data.frame(nimbleOut$summary$all.chains)

message('alpha = ',nimbleAC[grep('alpha',rownames(nimbleAC)),1])
message('beta = ',nimbleAC[grep('beta',rownames(nimbleAC)),1])
message('theta = ',nimbleAC[grep('theta',rownames(nimbleAC)),1])

x = (obsdata$timestamp_0-obsdata$timestamp_0[1])/3600


plot(x,nimbleData$y,col='blue',lty=1,ylab='DNA conc (copies/μL)',xlab='Time since start of sampling (hr)',log='y')
#points(obsdata$timestamp_1,nimbleAC[grep('^x\\[',rownames(nimbleAC)),1],col='red',lty=1)
lines(x,nimbleAC[grep('^x\\[',rownames(nimbleAC)),1],col='red',lty=1)
y = nimbleAC[grep('^x\\[',rownames(nimbleAC)),1]
sd = nimbleAC[grep('sigma\\[',rownames(nimbleAC)),1]
yplus = y+sd
ymin = y-sd
ymin[ymin<=0] <- 1
segments(x,ymin,x,yplus,col='red')

plot(obsdata$timestamp_1,nimbleData$y,col='blue',lty=1,ylab='DNA conc (copies/μL)',xlab='Timestamp')
points(obsdata$timestamp_1,nimbleAC[grep('x_s\\[',rownames(nimbleAC)),1],col='red',lty=1)
lines(nimbleConsts$Mod_ts,nimbleAC[grep('x\\[',rownames(nimbleAC)),1],col='red',lty=1)
segments(x, y-sd,x,y+sd,col='red')
