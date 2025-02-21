### example of how to compare estimates of population trends while accounting for the precision of those trend estimates

library(googledrive)
library(tidyverse)
library(jagsUI)


# Setup 
# Loading the BBS trends 1978-2022 for a collection of Prairie species

trends <- read_csv("data/prairie_sp_trends1978_2022.csv")

# ### selecting out the continental long-term trends and the trend-relevant columns
tr1 = trends %>%
  filter(region %in% c("23","22"),
         BirdName %in% c("American Crow",
                        "Henslow's Sparrow",
                        "Eastern Meadowlark",
                        "Ring-necked Pheasant",
                        "Eastern Bluebird",
                        "Sedge Wren",
                        "Tree Swallow",
                        "Grasshopper Sparrow",
                        "Common Yellowthroat",
                        "Field Sparrow",
                        "Blue-winged Teal",
                        "Savannah Sparrow",
                        "Northern Carinal",
                        "American Goldfinch",
                        "Clay-colored Sparrow",
                        "Bobolink",
                        "Killdeer",
                        "Eastern Kingbird",
                        "Red-winged Blackbird",
                        "Brown-headed Cowbird",
                        "Song Sparrow",
                        "Vesper Sparrow",
                        "Common Grackle",
                        "Barn Swallow",
                        "Mourning Dove",
                        "Western Meadowlark")) %>%
  select(BirdName,trend, trend_q_0.025, trend_q_0.975,n_routes,region) %>% 
  rename(species = BirdName,
         trend_lci = trend_q_0.025,
         trend_uci = trend_q_0.975)







### funciton to transform %-change trends to log-scale geometric mean rates of change
log_trans <- function(x){
  log((x/100)+1)
}

## transforming the 95% CIs into an estimate of precision (1/variance)
##pivoting wider to separate two BCR-specific trend estimates
##using only species that have estimates for both, of course
   trlong <- tr1 %>% mutate(betahat = log_trans(trend),
           log_lci = log_trans(trend_lci),
           log_uci = log_trans(trend_uci)) %>%
    mutate(varhat = ((log_uci-log_lci)/(1.96*2))^2) %>% 
     pivot_wider(id_cols = species,names_from = region, values_from = c(betahat,varhat,n_routes) ) %>% 
     filter(.,complete.cases(betahat_22,betahat_23))

 
   nspecies = nrow(trlong)
   nroutes = as.matrix(trlong[,c("n_routes_22","n_routes_23")])
   betahat = as.matrix(trlong[,c("betahat_22","betahat_23")])
   varhat = as.matrix(trlong[,c("varhat_22","varhat_23")])
   
   jags_data <- list(nspecies = nspecies,
                     nroutes = nroutes,
                     betahat = betahat,
                     varhat = varhat)
# model -------------------------------------------------------------------

 modl <- "

### model for comparing trend estimates from two (indexed by e) different analyses (or regions), for species (indexed by s)
### for this comparison, e==1 are estimates from BCR 22 and e==2 from BCR 23
### Model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751
### in essence, this is a paired t-test style comparison, that accounts for the imprecision in each s*e trend estimate.

### input data compiled in R, consist of nspecies (the number of species), as well as 3 matrices: varhat, betahat, and n, each of which has nspecies rows and 2 columns
## varhat = estimates of the log-scale variances of trends
## betahat = estimates of the log-scale trends
## nroutes = sample size of the trends - number of routes or survey-sites used to generate the trend)
## nspecies = number of species

model{

for (e in 1:2) {

for(s in 1:nspecies) {

	varhat[s,e] ~ dgamma(p[s,e],lam[s,e]) #varhat = data = variance of trend estimates
	p[s,e] <- nroutes[s,e] / 2
	lam[s,e] <- p[s,e] * tau.betahat[s,e]
	tau.betahat[s,e] ~ dgamma(0.001,0.001)
	sd.betahat[s,e] <- 1/sqrt(tau.betahat[s,e])
	### above is prior for the estimated species level precision (tau-betahat), accounting fo rhte number of routes included
	
	betahat[s,e] ~ dnorm(beta[s,e],tau.betahat[s,e]) #betahat = data = trend estimates
	beta[s,e] ~ dnorm(0.0,1)
	#### beta is the prior for the estimated species level trend, accounting for the precision
	
	 pos[s,e] <- step(beta[s,e]) #tracks whether the trend is positive
	 neg[s,e] <- 1-pos[s,e] #tracks opposite
	

} # end of first s-species loop
  numpos[e] <- sum(pos[,e]) # if pos is tracked, then this estimates the number of species with positive trends
  numneg[e] <- sum(neg[,e]) # if neg is tracked, same

	} #end of e loop (indexing two models being compared)

### species level differences between each trend-source
### if dif[s] is positive = first-source of trend is more positive than second-source of trend
for(s in 1:nspecies) {
dif[s] <- ((exp(beta[s,1])-1)*100)-((exp(beta[s,2])-1)*100) # dif is a vector of the species-specific trend differences (%/year) after accounting for the imprecision of each estimate's trend and the group (survey/monitoring program) structure
difvar[s] <- sd.betahat[s,1]-sd.betahat[s,2] # difvar is a vector of the species-specific differences in sd of the trends (log-scale, not %/year) (not necessarily interesting)
} # end of second s-species loop


m.dif <- mean(dif[])
pos_m.dif <- step(m.dif) #tracks whether the trend is positive (posterior mean of this value is the probability that m.dif is positive)
m.difvar <- mean(difvar[])

dif.numpos <- numpos[1]-numpos[2]
dif.numneg <- numneg[1]-numneg[2]


} # end of model

   "  
   
   trend_comp = "trend_comparison.txt"
cat(modl,file = trend_comp)   


params <- c("m.dif",
            "m.difvar",
            "pos_m.dif",
            "beta",
            "dif",
            "dif.numneg")


burnInSteps = 20000            # Number of steps to "burn-in" the samplers. this is sufficient for testing, but you'll want to increase this
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



summr = as.data.frame(out$summary) %>%  #table showing a summary of the posterior distribution and some basic convergence stats for all the monitored parameters
  rownames_to_column()

sp_names <- trlong %>% 
  select(species)
# species differences in trends between the two sets
species_differences <- summr %>% 
  filter(grepl("dif[",rowname,fixed = TRUE)) %>% 
  bind_cols(sp_names)
  
  
species_differences

# alternative hyperparameter model -------------------------------------------------------------------

modl <- "

### same as above, but assuming that each species trend shares some underlying process that's driving the survey-specific trends

model{

for (e in 1:2) {

for(s in 1:nspecies) {

	varhat[s,e] ~ dgamma(p[s,e],lam[s,e]) #varhat = data = variance of trend estimates
	p[s,e] <- nroutes[s,e] / 2
	lam[s,e] <- p[s,e] * tau.betahat[s,e]
	tau.betahat[s,e] ~ dgamma(0.001,0.001)
	sd.betahat[s,e] <- 1/sqrt(tau.betahat[s,e])
	### above is prior for the estimated species level precision (tau-betahat), accounting fo rhte number of routes included
	
	betahat[s,e] ~ dnorm(beta[s,e],tau.betahat[s,e]) #betahat = data = trend estimates
	beta[s,e] ~ dnorm(mu[e],tau.beta[e]) #centering each species estimate around the trend-source-specific hyperparameter mu[e], 
	#### beta is the prior for the estimated species level trend, accounting for the precision
	
	 pos[s,e] <- step(beta[s,e]) #tracks whether the trend is positive
	 neg[s,e] <- 1-pos[s,e] #tracks opposite
	

} # end of first s-species loop
  numpos[e] <- sum(pos[,e]) # if pos is tracked, then this estimates the number of species with positive trends
  numneg[e] <- sum(neg[,e]) # if neg is tracked, same


mu[e] ~ dnorm(0,1) #mean national trend hyperparameter
mu_ppy[e] <- ((exp(mu[e])-1)*100)
tau.beta[e] ~ dgamma(0.001,0.001) # precision
sd.beta[e] <- 1/sqrt(tau.beta[e]) # SD among species in their source-specific trend estimates


	} #end of e loop (indexing two models being compared)

### species level differences between each trend-source
### if dif[s] is positive = first-source of trend is more positive than second-source of trend
for(s in 1:nspecies) {
dif[s] <- ((exp(beta[s,1])-1)*100)-((exp(beta[s,2])-1)*100) # dif is a vector of the species-specific trend differences (%/year) after accounting for the imprecision of each estimate's trend and the group (survey/monitoring program) structure
difvar[s] <- sd.betahat[s,1]-sd.betahat[s,2] # difvar is a vector of the species-specific differences in sd of the trends (log-scale, not %/year) (not necessarily interesting)
} # end of second s-species loop

m.dif <- mu_ppy[1]-mu_ppy[2] #now represents the difference between the hyperparameters
pos_m.dif <- step(m.dif) #tracks whether the trend is positive (posterior mean of this value is the probability that m.dif is positive)

m.dif_s <- mean(dif[]) #mean of species level differences
m.difvar_s <- mean(difvar[])

dif.numpos <- numpos[1]-numpos[2]
dif.numneg <- numneg[1]-numneg[2]


} # end of model

   "  

trend_comp = "trend_comparison_hyperparameters.txt"
cat(modl,file = trend_comp)   


params <- c("m.dif",
            "mu_ppy",
            "pos_m.dif",
            "dif",
            "dif.numneg")


burnInSteps = 20000            # Number of steps to "burn-in" the samplers. this is sufficient for testing, but you'll want to increase this
nChains = 3                   # Number of chains to run.
numSavedSteps=1000         # Total number of steps in each chain to save. this is sufficient for testing, but you'll want to increase this
thinSteps=50                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.



out_hyp = jagsUI(data = jags_data,
             parameters.to.save = params,
             n.chains = 3,
             n.burnin = burnInSteps,
             n.thin = thinSteps,
             n.iter = nIter,
             parallel = T,
             model.file = trend_comp)



summr_hyp = as.data.frame(out_hyp$summary) %>%  #table showing a summary of the posterior distribution and some basic convergence stats for all the monitored parameters
  rownames_to_column()


# species differences in trends between the two sets
species_differences_hyp <- summr_hyp %>% 
  filter(grepl("dif[",rowname,fixed = TRUE)) %>% 
  bind_cols(sp_names)

# mean difference
mean_diffs <- summr_hyp %>% 
  filter(rowname %in% c("mu_ppy[1]","mu_ppy[2]","m.dif","pos_m.dif")) 

mean_diffs

