This directory provides the output files from the NA Sampler and NA Bayes runs used in the analysis for the 2024 paper. 

na.sum is a summary file for the NA Sampler run. 

na.nad is the direct access file created by the NA Sampler package and lists all of the forward models sampled and their associated misfits.

The na.nad file is the required input for the NA Bayes package, which calculates Bayesian statistics describing the parameter estimation.
  
nab.out is the output of the NA Bayes program and can be used to plot marginal and joint posterior probability distributions, 
posterior covariances, resolution, etc for the inversion. See Sambridge (1999a,b) for more information. 


References:
1. Sambridge, M. (1999a). Geophysical inversion with a neighbourhood algorithm—I. Searching a parameter space. Geophysical journal international, 138(2), 479-494.
2. Sambridge, M. (1999b). Geophysical inversion with a neighbourhood algorithm—II. Appraising the ensemble. Geophysical Journal International, 138(3), 727-746.
