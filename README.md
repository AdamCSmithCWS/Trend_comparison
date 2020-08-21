# Trend_comparison
example Bayesian hierarchical model for comparing species' population trends between two regions



Model for comparing trend estimates from two (indexed by e) different analyses (or regions), for species (indexed by s)
for this comparison, e==1 are estimates from the US and e==2 from Canada.
The model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751. In essence, this is a paired t-test style comparison, that accounts for the imprecision in each trend estimate.


## Adding a correlation-estimation model
example of a similar model for estimating the correlation between two collections of trend. For example, comparing regional trends from BBS and eBird, or comparing short-term and long-term trends for a single species from the BBS

example_script_correlation.R includes an example of this using Barn Swallow data for short and long-term BBS trends.

Model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751 and includes covariance components and corelation coefficient aspects from Matzke et al. 2017, https://doi.org/10.1525/collabra.78







