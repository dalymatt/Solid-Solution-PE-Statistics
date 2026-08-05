# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:32:07 2026

@author: abask
"""

import numpy as np
import pandas as pd
import copy
from collections import Counter
import rdf_coord as rc
import Potential as pot
from read_write_data import load_obj


#------Cut off Radius---------#
rcut=5.6 #FeNiCr
#rcut=5.8037500000 #FeNiCrCoCu


#------EAM/Alloy Potential---------#
fname='FeNiCr.eam.alloy'
#fname='FeNiCrCoCu-with-ZBL.eam.alloy'

#------Lattice Parameter---------#
#$$$$Random Systems$$$$$$
#ao=3.51036 #FeNiCr 73 8 19 Random
#ao=3.54939 # equimolar FeNiCrCoCu
#ao=3.4986 #Pure Fe in Fe Ni Cr
#ao=3.53073 # equimolar NiCrCo
#ao=3.546 # equimolar NiCrCoCu

#$$$$SRO Systems$$$$$$
#ao=3.50931 #FeNiCr 73 8 19 SRO alpha 0.05
ao=3.5150 #FeNiCr 73 8 19 SRO alpha -0.05


#------Compositions of the systems---------#
comp=np.array([0.73,0.08,0.19]) # FeNiCr
#comp=np.array([0.2,0.2,0.2,0.2,0.2]) # FeNiCrCoCu
#comp=np.array([0, 0.333,0.333,0.334, 0]) # NiCrCo
#comp=np.array([0, 0.25 ,0.25 ,0.25, 0.25]) # NiCrCoCu
#comp=np.array([1,0,0])

# 
fname_FCC = 'Fe_Perf.41.lammps'

#------W-C Parameters for upto 5 neighbor shell---------#
#alpha_file = 'alpha_FeNiCr_SS_point05.npy'
alpha_file = 'alpha_FeNiCr_SS_minus_point05.npy'

#------Ordering COnditions---------#
use_alpha = True   # True for SRO systems
#use_alpha = False   # False for Random systems



# =========================================================
# LAMMPS DUMP READER (MATT-STYLE)
# =========================================================
def Lammps_atoms(fname):
    atoms = pd.read_csv(fname, skiprows=8, delimiter=' ', engine='python')
    atoms = atoms.rename(columns={
        "ITEM:": "id", "ATOMS": "x", "id": "y",
        "x": "z", "y": "type",
        "z": "c_csym2", "type": "c_eng2"
    })
    if 'c_csym' in atoms.columns:
        atoms.drop(columns=['c_csym', 'c_eng'], inplace=True)
    atoms = atoms.rename(columns={"c_csym2": "c_csym", "c_eng2": "c_eng"})
    return atoms[['x', 'y', 'z']].values


# =========================================================
# 1) PERFECT FCC STRUCTURE FACTOR
# =========================================================
atoms_FCC = Lammps_atoms(fname_FCC)
min_spot = np.argmin(np.sum(atoms_FCC**2, axis=1))

_, cn_FCC = rc.rdf_coord(atoms_FCC, rcut, atoms_FCC[min_spot])
cn_FCC[:,0] = cn_FCC[:,0] / cn_FCC[1,0]

#[rdf_FCC, cn_FCC] = rc.rdf_coord_fcc(ao,rcut)


# =========================================================
# 2) READ POTENTIAL + FCC COHESIVE ENERGY
# =========================================================
rrange, rhorange, rho, Fr, Pp = pot.potential_read(fname)

alpha = np.load(alpha_file)

# --- safety checks ---
nshell = cn_FCC.shape[0]
nelem  = len(comp)

if use_alpha:   # <-- control switch
    alpha = np.load(alpha_file)

    assert alpha.shape == (nshell, nelem, nelem), \
        f"alpha shape mismatch: expected {(nshell,nelem,nelem)}, got {alpha.shape}"
else:
    alpha = np.zeros((nshell, nelem, nelem))




cn_perf = copy.deepcopy(cn_FCC)
cn_perf[:,0] = cn_perf[:,0] * ao

coh_FCC, _, E_element_FCC = pot.potential_stats(
    rrange, rhorange, rho, Fr, Pp, comp, cn_perf, alpha
)


# =========================================================
# 3) VFE ENVIRONMENTS — FROM Vac_rel.pkl (MATT-STYLE)
# =========================================================
import pickle

with open('cn_vac.pkl', 'rb') as f:
    cn_Vac1 = pickle.load(f)

cn_VFE_raw = cn_Vac1

cn_VFE = []
for cn_tmp in cn_VFE_raw:
    cn_norm = copy.deepcopy(cn_tmp)
    cn_norm[:, 0] = cn_norm[:, 0] / 3.49869654884664  # REQUIRED normalization
    cn_VFE.append(cn_norm)

counts_VFE = Counter([tuple(map(tuple, c)) for c in cn_VFE])
cn_VFE_unique = [np.array(k) for k in counts_VFE.keys()]
cn_VFE_freq   = list(counts_VFE.values())



# =========================================================
# 4) VFE — MATT-STYLE (LINE-BY-LINE)
# =========================================================
stats_VFE_atoms_analytical = np.zeros((len(cn_VFE_unique), 6))

E_VFE_Bar = 0.0
res1 = 0.0
res2 = 0.0

for i in np.arange(0,len(cn_VFE_unique)):
    cn=copy.deepcopy(cn_VFE_unique[i])
    cn[:,0]=cn[:,0]*ao
    [coh_site,_,E_element_VFE]=pot.potential_stats(rrange,rhorange,rho,Fr,Pp,comp,cn, alpha)    
    coh_site['E']=coh_site['E'].fillna(0)
    stats_VFE_atoms_analytical[i,0]=i
    stats_VFE_atoms_analytical[i,1]=coh_FCC['E'].iloc[0]
    stats_VFE_atoms_analytical[i,2]=coh_FCC['E'].iloc[1]
    stats_VFE_atoms_analytical[i,3]=coh_site['E'].iloc[0]
    stats_VFE_atoms_analytical[i,4]=coh_site['E'].iloc[1]
    covar= np.sum(comp*E_element_VFE*E_element_FCC)-(coh_site['E'].iloc[0]*coh_FCC['E'].iloc[0])
    stats_VFE_atoms_analytical[i,5]=covar
    E_VFE_Bar=E_VFE_Bar+(coh_site['E'].iloc[0]-coh_FCC['E'].iloc[0])*cn_VFE_freq[i]
    Sig_excess=coh_site['E'].iloc[1]*coh_site['E'].iloc[1] + coh_FCC['E'].iloc[1]*coh_FCC['E'].iloc[1]-2*coh_site['E'].iloc[1]*coh_FCC['E'].iloc[1]
    res1=res1+ cn_VFE_freq[i]*(Sig_excess)
    res2=res2+2*Sig_excess*(cn_VFE_freq[i]*(cn_VFE_freq[i]-1))/2     
     
stats_VFE_atoms_analytical=pd.DataFrame(stats_VFE_atoms_analytical, columns=['VFE atom site', 'Efcc average', 'Efcc stdev', 'Evfe average', 'Evfe stdev', 'covar'])

Sig_VFE=np.sqrt(res1+res2)

Sig_VFE_check=np.sum(np.array(cn_VFE_freq)*np.array(cn_VFE_freq)*\
                     (stats_VFE_atoms_analytical['Efcc stdev'].values*stats_VFE_atoms_analytical['Efcc stdev'].values+\
                         stats_VFE_atoms_analytical['Evfe stdev'].values*stats_VFE_atoms_analytical['Evfe stdev'].values-\
                             2*stats_VFE_atoms_analytical['Evfe stdev'].values*stats_VFE_atoms_analytical['Efcc stdev'].values))

Sig_VFE_check=np.sqrt(Sig_VFE_check)


stats_VFE_atoms_analytical = pd.DataFrame(
    stats_VFE_atoms_analytical,
    columns=[
        'VFE atom site',
        'Efcc average',
        'Efcc stdev',
        'Evfe average',
        'Evfe stdev',
        'covar'
    ]
)


# =========================================================
# 5) TS ENVIRONMENTS
# =========================================================
import pickle

with open('TS_rel.pkl', 'rb') as f:
    cn_Mig1 = pickle.load(f)

cn_Mig_raw = cn_Mig1

cn_TS = []
for cn_tmp in cn_Mig_raw:
    cn_norm = copy.deepcopy(cn_tmp)
    cn_norm[:, 0] = cn_norm[:, 0] / 3.49869654884664   # <-- IDENTICAL normalization as VFE
    cn_TS.append(cn_norm)

counts_TS = Counter([tuple(map(tuple, c)) for c in cn_TS])
cn_TS_unique = [np.array(k) for k in counts_TS.keys()]
cn_TS_freq = list(counts_TS.values())



# =========================================================
# 6) TS — MATT-STYLE (LINE-BY-LINE)
# =========================================================
stats_TS_atoms_analytical = np.zeros((len(cn_TS_unique), 6))

E_TS_Bar = 0.0
E_VME_bar = 0.0
res1 = 0.0
res2 = 0.0

for i in range(len(cn_TS_unique)):

    cn_loc = copy.deepcopy(cn_TS_unique[i])
    cn_loc[:,0] = cn_loc[:,0]*ao

    coh_site, _, _ = pot.potential_stats(
        rrange, rhorange, rho, Fr, Pp, comp, cn_loc, alpha
    )

    stats_TS_atoms_analytical[i,0] = i
    stats_TS_atoms_analytical[i,1] = coh_FCC['E'].iloc[0]
    stats_TS_atoms_analytical[i,2] = coh_FCC['E'].iloc[1]
    stats_TS_atoms_analytical[i,3] = coh_site['E'].iloc[0]
    stats_TS_atoms_analytical[i,4] = coh_site['E'].iloc[1]
    stats_TS_atoms_analytical[i,5] = 0.0

    Ni = cn_TS_freq[i]

    E_TS_Bar = (
        E_TS_Bar
        + (coh_site['E'].iloc[0] - coh_FCC['E'].iloc[0]) * Ni
    )
    
    E_VME_bar = E_TS_Bar - E_VFE_Bar

    Sig_excess = (
        coh_site['E'].iloc[1]*coh_site['E'].iloc[1]
        + coh_FCC['E'].iloc[1]*coh_FCC['E'].iloc[1]
        - 2*coh_site['E'].iloc[1]*coh_FCC['E'].iloc[1]
    )

    res1 = res1 + Ni * Sig_excess
    res2 = res2 + 2 * Sig_excess * (Ni*(Ni-1)) / 2

Sig_TS = np.sqrt(res1 + res2)
Sig_VME = np.sqrt(Sig_TS**2 + Sig_VFE**2 )

stats_TS_atoms_analytical = pd.DataFrame(
    stats_TS_atoms_analytical,
    columns=[
        'TS atom site',
        'Efcc average',
        'Efcc stdev',
        'ETS average',
        'ETS stdev',
        'covar'
    ]
)

print("E_VFE  =", E_VFE_Bar)
print("Sig_VFE=", Sig_VFE)
print("E_VME  =", E_VME_bar)
print("Sig_VME=", Sig_VME)
