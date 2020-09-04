library(bbsBayes)
library(tidyverse)
library(plyr)
library(jagsUI)
library(ggmcmc)
library(beepr)

#set up for WOTH

#load eBird simulated TRUE trends
sim_trends <- read.csv("data/ebird-trends_simulated.csv", skipNul = TRUE,colClasses = c("character","integer",rep("numeric",9)))

#will need to load eBird estimated trends as well, should use same methods as the sim for formatting
#not needed for model selection


#functions
#log transform fxn trends and CI-values  
log_trans <- function(x){
  log((x/100)+1)
}

# function to transform CI width into variance ----------------------------
## this function will work with trends on the %-change or log-scale
CI_to_var <- function(lci, #value of lower credible limit
                      uci, #value of upper credible limit
                      CI_p){ #95% CI in data, if otherwise, insert the type of CI here (e.g., if 80% CI, CI_p = 0.8)
  qth <- qnorm(1-((1-CI_p)/2))
  return(((uci-lci)/(qth*2))^2) #returns the variance
}


for(h in 1:10){

#load the simulated eBird data run through our models
wothsim <- read.csv(paste("data/eBird_WOTH_sim", h, "_spatial_CAR_trends.csv", sep = ""))



#put TRUE trends and MOD trends into one df

#all steps combined - do this for each dataset, then stack dataframes and pivot

#WOTH scenario 1 true trends
#use same code for eBirdTM estimates
sim_trends_true <- sim_trends %>% 
  mutate(Region = paste(trunc(lat), trunc(lon), sep = "_")) %>%   #stratify to region
  mutate(betahat = log_trans(abd_ppy_true),       #log transform trends and CIs
         log_lci = log_trans(abd_ppy_lower),
         log_uci = log_trans(abd_ppy_upper)) %>%
  mutate(tau.betahat = 1/CI_to_var(lci = log_lci,uci = log_uci, CI_p = 0.8)) %>%  #calculate precision
  filter(rss_name == "woothr_breeding",                #filter for species and scenario
         scenario == 1) %>%  ################  this is where we could make it a function for all scenarios
  group_by(Region) %>%
  dplyr::summarize(betahat = mean(betahat),
                   tau.betahat = mean(tau.betahat)) #
#%>% #summarize by region (stratum)
sim_trends_true <- mutate(sim_trends_true,dataset = rep("TRU", length(sim_trends_true$Region)))   #add column with dataset identification
#474 obs

  
#WOTH scenario 1 estimated trends using our model and ebird data
sim_trends_mod <- wothsim %>% 
  mutate(dataset = rep("MOD", length(wothsim$Region))) %>% #add column with dataset identification
  mutate(betahat = log_trans(Trend),        # log transform
         log_lci = log_trans(Trend_Q0.025),
         log_uci = log_trans(Trend_Q0.975)) %>%
  mutate(tau.betahat = 1/CI_to_var(lci = log_lci,uci = log_uci, CI_p = 0.95)) # precision
#436 obs


#stack dataframes, create full dataframe for scenario one with both true trends and estimated trends

sc1 <- rbind.fill(sim_trends_true, sim_trends_mod) %>%
  pivot_wider(id_cols = Region, names_from = dataset, values_from = c(betahat,tau.betahat) ) %>% 
  filter(.,complete.cases(betahat_TRU, betahat_MOD))

#note: for final summaries, will need to transform back from the log scale, but not necessary for model selection



nstrata = nrow(sc1)
betahat = as.matrix(sc1[,c("betahat_TRU","betahat_MOD")])
tau.betahat = as.matrix(sc1[,c("tau.betahat_TRU","tau.betahat_MOD")])

jags_data <- list(nstrata = nstrata,
                  betahat = betahat,
                  tau.betahat = tau.betahat)

M1 = lm(data = sc1,
        formula = betahat_TRU~betahat_MOD,
        weights = tau.betahat_MOD)

### model for comparing trend estimates from two (indexed by e) different analyses, for different regions or strata (indexed by s)
### for this comparison, e==1 are estimates from the BBS and e==2 from eBird
### Model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751
### and includes covariance components and corelation coefficient aspects from Matzke et al. 2017, https://doi.org/10.1525/collabra.78

### in essence, this is a paired t-test style comparison, that accounts for the imprecision in each s*e trend estimate.

### input data compiled in R, consist of nstrata (the number of strata), as well as 2 matrices: tau.betahat and betahat; 
###both have nstrata rows and 2 columns
## tau.betahat = estaimtes of the precision of the trend estimates (1/variance)
## betahat = estimates of the trends
## nstrata = number of strata
modl <- "

model{
  
  for (e in 1:2) {
    
    for(s in 1:nstrata) {
      
      betahat[s,e] ~ dnorm(beta[s,e],tau.betahat[s,e]) #betahat = data = trend estimates, tau.betahat = data = precision of trend estimate
      
    } # end of s loop
    
    mu[e] ~ dnorm(0,1) #mean trend for each survey
  sd.beta[e] <- 1/sqrt(tau.beta[e])
    
    
  } #end of e loop (indexing two models being compared)
      tau.beta[1] ~ dscaled.gamma(0.001,50) 
  tau.beta[2] ~ dscaled.gamma(0.01,50) 
  
  
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
"rho", #rho is the correlation coefficient
"mu")

burnInSteps = 5000        #2000    # Number of steps to "burn-in" the samplers. this is sufficient for testing, but you'll want to increase this
nChains = 3                   # Number of chains to run.
numSavedSteps=2000   #1000      # Total number of steps in each chain to save. this is sufficient for testing, but you'll want to increase this
thinSteps=10          #10         # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.



start.time <- Sys.time()

out = jagsUI(data = jags_data,
parameters.to.save = params,
n.chains = 3,
n.burnin = burnInSteps,
n.thin = thinSteps,
n.iter = nIter,
parallel = T,
model.file = trend_comp)
beep("mario")

end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken


summr = out$summary #table showing a summary of the posterior distribution and some basic convergence stats for all the monitored parameters
 
out_ggs = ggs(out$samples)
#ggmcmc(out_ggs,file = "correlation_convergence_summaries.pdf", param_page = 8)
ggmcmc(out_ggs,file = paste("correlation_convergence_summaries_WOTH", h, ".pdf", sep = ""), param_page = 8)



cor <- cor(jags_data$betahat[,1],jags_data$betahat[,2]) # raw Pearson correlation without considering uncertainty
rho <- out$summary["rho",] ## estimated correlation accounting for uncertainty
cor
rho

save(list = c("cor", "rho"), file = paste("WOTH_sim", h, "results.RData", sep = ""))


  
