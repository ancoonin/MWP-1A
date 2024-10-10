This directory provides the output files from the NA Sampler and NA Bayes runs used in the analysis for the 2024 paper. 

na.sum is a summary file for the NA Sampler run. 

na.nad is the direct access file created by the NA Sampler package and lists all of the forward models sampled and their associated misfits.

The na.nad file is the required input for the NA Bayes package, which calculates Bayesian statistics describing the parameter estimation.
  
nab.out is the output of the NA Bayes program and can be used to plot marginal and joint posterior probability distributions, 
posterior covariances, resolution, etc for the inversion. See Sambridge (1999a,b) for more information. 

prem.l75.umvm7.lmvm7_elastic contains the outputs of the MWP-1A inversion assuming an elastic only solid Earth response
prem.l75.umvm7.lmvm7_maxwell contains the outputs of the MWP-1A inversion assuming a Maxwell Earth rheology
prem.l75.umvm7.lmvm7_EXB contains the outputs of the MWP-1A inversion assuming an extended burgers solid Earth rheology
prem.l75.umvm7.lmvm7_mc2001 contains the outputs of the MWP-1A inversion, assuming the solid Earth deformation to follow
the empirical transient rheology law of McCarthy et al. (2011)

All of the viscoelastic Earth models assume the background structure of the VM7 Earth model (Roy and Peltier, 2018).

References:
1. Sambridge, M. (1999a). Geophysical inversion with a neighbourhood algorithm—I. Searching a parameter space. Geophysical journal international, 138(2), 479-494.
2. Sambridge, M. (1999b). Geophysical inversion with a neighbourhood algorithm—II. Appraising the ensemble. Geophysical Journal International, 138(3), 727-746.
3. McCarthy, C., Takei, Y., & Hiraga, T. (2011). Experimental study of attenuation and dispersion over a broad frequency range: 2. The universal scaling of polycrystalline materials. Journal of Geophysical Research: Solid Earth, 116(B9).
4. Roy, K., & Peltier, W. R. (2018). Relative sea level in the Western Mediterranean basin: A regional test of the ICE-7G_NA (VM7) model and a constraint on late Holocene Antarctic deglaciation. Quaternary Science Reviews, 183, 76-87.
