#!/usr/bin/env python 

import numpy
from numpy import array
import dadi
import sys
import pylab 
import Demographics2D

"""
Script for optimization of model fit to dadi input files. For each tidal marsh-interior comparison
need to change the input file for Misc.make_data_dict, population sizes, and output file. 

This will run optimization 10 times. A final optimization is run with the parameters from the 
optimization with the highest log-likelihood 
"""

#parse snps file
dd = dadi.Misc.make_data_dict("AnuInt.final.dadi")
fs = dadi.Spectrum.from_data_dict(dd, ['INT', 'ANU'], [16,13],
                                       polarized=False)
                                       
print fs.S()
                                

# model and data settings                                      
data = fs
ns = data.sample_sizes
pts_l = [40,50,60]
                                      
#for each model I will have different levels of:
#model name

nu1 = 5
nu2 = 5
Tsp = 1
m1 = 0.1
m2 = 0.1


params = [nu1, nu2, Tsp, m1,m2]

upper_bounds = [100, 500, 100, 100, 100]

lower_bounds = [1e-2, 1e-2, 1e-3, 0,0]


IMmod = Demographics2D.split_mig
func = IMmod
func_ex = dadi.Numerics.make_extrap_log_func(func)
maxiters=10
OutFileName = "./AnulusInt.txt"
OutFile = open(OutFileName, 'a')   
for j in range(10):
	print(j)
	# Perturb  parameters before optimization. This does so by taking each
	# parameter a up to a factor of two up or down.
	p0 = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bounds,
								  lower_bound=lower_bounds)
		  
	# Do the optimization. By default we assume that theta is a free parameter,
	# since it's trivial to find given the other parameters. If you want to fix
	# theta, add a multinom=False to the call.
	print('Beginning optimization ************************************************')
	popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
									   lower_bound=lower_bounds,
									   upper_bound=upper_bounds,
									   verbose=len(p0), maxiter=maxiters)
	# The verbose argument controls how often progress of the optimizer should be
	# printed. It's useful to keep track of optimization process.
	print('Finshed optimization **************************************************')

	model = func_ex(popt, ns, pts_l)
	# Likelihood of the data given the model AFS.
	ll_opt = dadi.Inference.ll_multinom(model, data)
	theta = dadi.Inference.optimal_sfs_scaling(model,data)
	param_out = list(popt) + [ll_opt] + [theta]
	print(param_out)
	writeparam = str(param_out) + "\n"
	OutFile.write(writeparam)
OutFile.close()		