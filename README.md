# Trend_comparison
example Bayesian hierarchical model for comparing species' population trends between two regions or programs



Model for comparing trend estimates from two (indexed by e) different analyses (or regions), for species (indexed by s)
for this comparison, comparing BBS estimates from two regions (BCRs). Easily modified to compare estimates from one BCR with estimates from a local monitoring program within that BCR.

The model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751. In essence, this is a paired t-test style comparison, that accounts for the imprecision in each trend estimate.


## Alternate correlation-estimation model
In the alternates folder there is another example model that may be worth considering. It's a similar model for estimating the dfference between two sets of trends, while also accounting for the correlation between two collections of trend. For example, comparing regional trends from BBS and eBird, or comparing short-term and long-term trends for a single species from the BBS.
This model also uses estimates of the trend precision as data, instead of using the chi-squared distribution and the number of routes to estimate the precision.

example_script_correlation.R includes an example of this using Barn Swallow data for short and long-term BBS trends.

Model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751 and includes covariance components and corelation coefficient aspects from Matzke et al. 2017, https://doi.org/10.1525/collabra.78







