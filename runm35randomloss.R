

library(nimble)

source("./Models/m35randomloss.R")

moddata <- read.csv("./Data/dolphin_tide_minutefreq.csv")
obsdata <- read.csv("./Data/ESP_timestamps_mLseawater.csv")

nimbleData <- list(
  y = obsdata$copies_per_mLseawater
  
)

nimbleConsts <- list(
  nTSo = nrow(obsdata)
)

nimbleInits <- list(alpha = rep(0.5, nrow(obsdata)),
                    beta = 1,
                    xinit = 100,
                    x = rep(1, nrow(obsdata)),
                    sigma = rep(1, nrow(obsdata)),
                    theta = 1.0
)

nimbleParams <- list("x", "alpha", "beta", "xinit", "sigma","sd_log_x","theta")

model <- nimbleModel(code = m35randomloss,
                     data = nimbleData,
                     constants = nimbleConsts,
                     inits = nimbleInits,
                     check = FALSE)

set.seed(20240401) 

nimbleOut <- nimbleMCMC(model, 
                        monitors = nimbleParams, 
                        thin = 100, niter = 1000000, nburnin = 750000, nchains = 4,
                        summary = TRUE,WAIC = TRUE)

save(nimbleOut, file = "./Results/nimbleOut_m35randomloss.RData")

load("./Results/nimbleOut_m35randomloss.RData")

nimbleAC <- as.data.frame(nimbleOut$summary$all.chains)

message('theta = ',nimbleAC['theta',1])
#message('alpha = ',nimbleAC['alpha',1])
message('beta = ',nimbleAC['beta',1])
#message('gamma = ',nimbleAC['gamma',1])


x = (obsdata$timestamp_0-obsdata$timestamp_0[1])/3600


plot(x,nimbleData$y,col='blue',lty=1,ylab='DNA conc (copies/mL)',xlab='Time since start of sampling (hr)',log='y')
#points(obsdata$timestamp_1,nimbleAC[grep('^x\\[',rownames(nimbleAC)),1],col='red',lty=1)
lines(x,nimbleAC[grep('^x\\[',rownames(nimbleAC)),1],col='red',lty=1)
y = nimbleAC[grep('^x\\[',rownames(nimbleAC)),1]
sd = nimbleAC[grep('sigma\\[',rownames(nimbleAC)),1]
yplus = y+sd
ymin = y-sd
ymin[ymin<=0] <- 1
segments(x,ymin,x,yplus,col='red')

#plot(obsdata$timestamp_1,nimbleData$y,col='blue',lty=1,ylab='DNA conc (copies/μL)',xlab='Timestamp')
#points(obsdata$timestamp_1,nimbleAC[grep('x_s\\[',rownames(nimbleAC)),1],col='red',lty=1)
#lines(nimbleConsts$Mod_ts,nimbleAC[grep('x\\[',rownames(nimbleAC)),1],col='red',lty=1)
#segments(x, y-sd,x,y+sd,col='red')

df <- data.frame (
  hours_since_ESP_start = x,
  ESP_DNA_conc = nimbleData$y,
  model = nimbleAC[grep('^x\\[',rownames(nimbleAC)),1],
  sigma = sd
)
write.csv(df, "./Results/m35randomloss_plotting.csv")

