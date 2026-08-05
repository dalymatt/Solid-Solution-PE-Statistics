"""
Advanced Materials and Microstructures Lab 
University of Illinois at Chicago, 2022
https://amml.lab.uic.edu/

Please cite: Jagatramka et al. (2022), Comp. Mater. Sci. 214 111763, doi: 10.1016/j.commatsci.2022.111763

"""
#Here, necessary Python packages will be imported and other written scripts will be called to perform the energy statistics calculations herein
import numpy as np  #Python package used commonly for data manipulation involving arrays, matrices, etc
import Potential as pot #Performs potential energy calculations
import rdf_coord as rc #Generated coordination number (CN) and radial distribution function (RDF) relations based on lattice parameter and cutoff distance
#%%
#FeNiCr
# ao = 3.5225 #Lattice parameter value in Angstroms
# rcut = 5.6 #Cutoff distance (Angstroms) which will restrict coordination shells generated
#NiCo
ao=3.512 #Lattice parameter value in Angstroms
rcut=6.5 #Cutoff distance (Angstroms) which will restrict coordination shells generated
#%%
#Refers to the composition of the system,
#the dimensionality defines the number of components
comp=np.array([0.5,0.5])
# comp=np.array([0.73307849132,0.07540562983,0.19151587883]) #SS Composition for FeNiCr
# comp = np.array([0.33, 0.33, 0.34])  
#%%
#Generates RDF and CN for perfect lattice based on ao and rc values, 
#user can change fcc part to hcp or bcc based on material's lattice
[rdf, cn] = rc.rdf_coord_fcc(ao,rcut) 
#[rdf, cn] = rc.rdf_coord_hcp(ao,rcut) 
#%%
#EAM potential file, in setfl format 
fname='NiCo-lammps-2014.alloy'
# fname='FeNiCr.eam.alloy'
#%%
#Reading the EAM potential dataset declared above and extracting needed parameters from file
[rrange, rhorange, rho, Fr, Pp] = pot.potential_read(fname)
#%%
#Performs Cohesive Energy Calculations and Statistics based on chosen EAM file of chosen system
# E_element is solute-level peratom cohesive energy used for GPFE calculations see Eq (11) in the main text for further details
# form_E gives the Avg and Mean of the EAM terms and covar is the covariance for the final cohesive energy calculation
[form_E, covar,E_element] = pot.potential_stats(rrange, rhorange, rho, Fr, Pp, comp, cn)
print(form_E) # The units here are eV/atom 
#%%
#Generates RDF and CN for different faulted strucutre based on ao, rcut, and cn_FCC values, 
for ii in ['USF','ISF','UTF1','ESF','UTF2','TF']:
    [ globals()['rdf_%s' % ii], globals()['cn_%s' % ii]] = rc.rdf_coord_fault(ao,rcut,cn,ii) 
#%%
# Generates the final mean and standard deviation to the interplanar fault energies of the GPFE landscape
# The energy units are eV/atom
for ii in ['USF','ISF','UTF1','ESF','UTF2','TF']: 
    [ globals()['form_E_%s' % ii], globals()['covar_%s' % ii]] = pot.potential_stats_fault(rrange, rhorange, rho, Fr, Pp, comp, cn,globals()['cn_'+ii],form_E,E_element,ii)
    print (globals()['form_E_'+ii])  # The units here are eV/atom 
