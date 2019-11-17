#!/usr/bin/env python 

import numpy
from numpy import array
import dadi
import sys
import pylab 
import Demographics2D

Usage="""Script to use Godambe information matrix methods (Coffmann et al. 2016) to
estimate parameter uncertainties based on bootstrap sfs generated 
using script: SAVSphylogeog_bootstrap.py
"""

#parse snps file
dd = dadi.Misc.make_data_dict(""~/DadiAnalyses/BenhamCheviron_dadi_input_files/AnuInt.final.dadi"")
SAVSfoldedSpec = dadi.Spectrum.from_data_dict(dd, ['INT','ANU'], [16,13],
                                      polarized=False)
                                      
# model and data settings                                      
data = SAVSfoldedSpec
ns = data.sample_sizes
pts_l = [30,40,60]

#Optmized parameters for each population
Atratus_popt = [1.5658714646843674, 33.321172251082515, 13.130210559190406, 0.0027724026796928634, 0.36733285192332554]
Belding_popt = [5.1679212106714152, 37.660892758847581, 6.2542132991674952, 0.14252249296043207, 0.11866353696550071]
anulus_popt = [4.3756141453448274, 45.507960753977507, 7.2594896020183786, 0.081303530839243507, 0.097776299962358301]
Guttatus_popt = [6.4821771589202077, 78.188374211046707, 15.649773712862521, 2.1884006476798321e-05, 0.075682318044207358]
Magdalanae_popt = [15.756676509894548, 0.77788111604571453, 11.670525615204594, 0.80015922851904309, 0.00011262125789092772]
Humboldt_popt = [6.9878838323055215, 31.623349159151473, 5.5224880163472214, 1.9520996234098589, 0.84199046786216747]
GIWA_popt = [6.4572714561978239, 49.106352281790485, 9.0128353962974934, 1.928492691564363, 0.024205358232868197]
NAPA_popt = [26.444047584110983, 5.1676575592728353, 4.842020526168537, 0.25813553952495916, 2.1945290590054505]
SFB_popt = [46.91565477827038, 17.538289977368979, 10.152322900232317, 0.26512796732542621, 0.33624148355304839]
Morro_popt = [22.669182500333715, 3.3135380324428558, 3.7608227621715704, 0.25331399820204492, 0.25674119893339625]

IMconstant = Demographics2D.split_mig

OutFile = open('./CIanalyses/Uncerts_output_AnuInt.txt', 'w')  

func = IMconstant
func_ex = dadi.Numerics.make_extrap_log_func(func)

#requires bootstraps generated from BenhamCheviron_bootstrap.py as input in all_boot
all_boot = [dadi.Spectrum.from_file('./CIanalyses/bootstraps_anu/fsboot_{0:02d}'.format(ii)) 
			for ii in range(100)]
uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, SFB_popt, data, multinom=True, return_GIM=True)                       
OutFile.write('{0} estimated sd GIM uncerts: {1}'.format("IMconstant_anu",uncerts[0]))                            
inverseMat = numpy.linalg.inv(uncerts[1])

OutFile.close()