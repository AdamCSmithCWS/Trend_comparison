### example of how to estimate the correlation of a series of estimates of population trends while accounting for the precision of those trend estimates

library(googledrive)
library(tidyverse)
library(jagsUI)
library(ggmcmc)

############ an example set of trends from the 2018 BBS analysis
# Setup downloading the 2018 BBS trends from my Google Drive --------------
### this commented out section doesn't need to be re-run
# tdwn <- drive_download(file = "All 2018 BBS trends.csv",
#                        path = "data/All_2018_BBS_trends.csv")
# 
# 
# trends <- read.csv("data/All_2018_BBS_trends.csv",stringsAsFactors = F)


# ### selecting out the stratum-level long-term and short-term trends and the trend-relevant columns
tr1 = trends %>%
  filter(Region_type == "stratum", species == "Barn Swallow") %>%
  select(Trend_Time,Trend, Trend_Q0.025, Trend_Q0.975,Number_of_Routes,Region)

tr1$trend_type = gsub(pattern = "-term",replacement = "",x = tr1$Trend_Time)
write.csv(tr1,"data/Barn_Swallow_long_and_short_trends.csv")

### this tr1 dataset provides a simple set of trend estimates to be correlated
### the model below then estimates the correlation between the short and long-term trends in strata for a single specie




# reading in saved BARS BBS trends ---------------------------------------

tr = read.csv("data/Barn_Swallow_long_and_short_trends.csv")

### funciton to transform %-change trends to log-scale geometric mean rates of change
log_trans <- function(x){
  log((x/100)+1)
}

CI_p <- 0.95 #95% CI in data, if otherwise, insert the type of CI here (e.g., if 80% CI, CI_p <- 0.8)
qth <- qnorm(1-((1-CI_p)/2))
## transforming the 95% CIs into an estimate of precision (1/variance)
##pivoting wider to separate US and Canadian estimates
##filtering to species that have estimates for both
   trlong <- tr %>% mutate(betahat = log_trans(Trend),
           log_lci = log_trans(Trend_Q0.025),
           log_uci = log_trans(Trend_Q0.975)) %>%
    mutate(tau.betahat = 1/(((log_uci-log_lci)/(qth*2))^2)) %>% 
     pivot_wider(id_cols = Region,names_from = trend_type, values_from = c(betahat,tau.betahat) ) %>% 
     filter(.,complete.cases(betahat_Long,betahat_Short))

 
   nstrata = nrow(trlong)
   betahat = as.matrix(trlong[,c("betahat_Long","betahat_Short")])
   tau.betahat = as.matrix(trlong[,c("tau.betahat_Long","tau.betahat_Short")])
   
   jags_data <- list(nstrata = nstrata,
                     betahat = betahat,
                     tau.betahat = tau.betahat)
# model -------------------------------------------------------------------

 modl <- "

### model for comparing trend estimates from two (indexed by e) different analyses, for different regions or strata (indexed by s)
### for this comparison, e==1 are estimates from the BBS and e==2 from eBird
### Model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751
### and includes covariance components and corelation coefficient aspects from Matzke et al. 2017, https://doi.org/10.1525/collabra.78

### in essence, this is a paired t-test style comparison, that accounts for the imprecision in each s*e trend estimate.

### input data compiled in R, consist of nstrata (the number of strata), as well as 2 matrices: tau.betahat and betahat; both have nstrata rows and 2 columns
## tau.betahat = estaimtes of the precision of the trend estimates (1/variance)
## betahat = estimates of the trends
## nstrata = number of strata

model{

for (e in 1:2) {

for(s in 1:nstrata) {

	betahat[s,e] ~ dnorm(beta[s,e],tau.betahat[s,e]) #betahat = data = trend estimates, tau.betahat = data = precision of trend estimate

} # end of s loop

mu[e] ~ dnorm(0,1) #mean trend for each survey
tau.beta[e] ~ dscaled.gamma(1,10) 
sd.beta[e] <- 1/sqrt(tau.beta[e])


	} #end of e loop (indexing two models being compared)


#### multivariate normal structure for the among-strata between survey variation in trends

for(s in 1:nstrata) {

beta[s,1:2] ~ dmnorm(mu[1:2],ISigma_cov[1:2,1:2])
}

rho ~ dunif(-1,1) #estimated Pearson correlation coefficient

### sigma_cov is the covariance matrix
Sigma_cov[1,1] <- pow(sd.beta[1],2)
Sigma_cov[1,2] <- rho*sd.beta[1]*sd.beta[2]
Sigma_cov[2,1] <- rho*sd.beta[1]*sd.beta[2]
Sigma_cov[2,2] <- pow(sd.beta[2],2)
ISigma_cov[1:2,1:2] <- inverse(Sigma_cov[1:2,1:2])


for(s in 1:nstrata) {
dif[s] <- beta[s,1]-beta[s,2] # dif is a vector of the strata-specific trend differences after accounting for the imprecision of each estimate's trend and the group (survey/monitoring program) structure
} # end of second s-strata loop


m.dif <- mean(dif[])

} # end of model

   "  
   
   trend_comp = "trend_comparison_correlation.txt"
cat(modl,file = trend_comp)   


params <- c("m.dif", #this is the mean difference (trend1 - trend2) it's an estimate of the difference in the average fine-scale trends (like a difference in teh continental estimate accounting for uncertainty of the local trends, but ignoring the abundance-weights)
            "Sigma_cov",#covariance matrix
            #"dif", # these values could be mapped to show the spatial pattern in the differences between trends
            #"beta", #these beta values could be used as posterior estimates of the regional trends aftern accounting for the precision and correlation
            "rho") #rho is the correlation coefficient


burnInSteps = 2000            # Number of steps to "burn-in" the samplers. this is sufficient for testing, but you'll want to increase this
nChains = 3                   # Number of chains to run.
numSavedSteps=1000         # Total number of steps in each chain to save. this is sufficient for testing, but you'll want to increase this
thinSteps=10                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.



out = jagsUI(data = jags_data,
             parameters.to.save = params,
             n.chains = 3,
             n.burnin = burnInSteps,
             n.thin = thinSteps,
             n.iter = nIter,
             parallel = T,
             model.file = trend_comp)



summr = out$summary #table showing a summary of the posterior distribution and some basic convergence stats for all the monitored parameters

out_ggs = ggs(out$samples)
ggmcmc(out_ggs,file = "correlation_convergence_summaries.pdf", param_page = 8)



cor(jags_data$betahat[,1],jags_data$betahat[,2]) # raw Pearson correlation without considering uncertainty
out$summary["rho",] ## estimated correlation accounting for uncertainty


