# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:49:52 2026

@author: abask
"""

import numpy as np

from vfe_vme_functions import (
    print_summary,
    run_vfe_vme_calculation,
)


# =========================================================
# SYSTEM INPUTS
# =========================================================

POTENTIAL_FILE = "FeNiCr.eam.alloy"

LATTICE_PARAMETER = 3.51036

COMPOSITION = np.array([
    0.73,
    0.08,
    0.19,
])


# =========================================================
# PRE-EXPORTED COORDINATION FILES
# =========================================================

FCC_COORDINATION_FILE = "cn_FCC.pkl"

VACANCY_ENVIRONMENT_FILE = "cn_vac.pkl"

TRANSITION_STATE_ENVIRONMENT_FILE = "TS_rel.pkl"


# =========================================================
# ORDERING CONDITION
# =========================================================

USE_ALPHA = True

ALPHA_FILE = "alpha_FeNiCr_SS_minus_point05.npy"

# For a random alloy:
# USE_ALPHA = False
# ALPHA_FILE = None


# =========================================================
# NORMALIZATION
# =========================================================

ENVIRONMENT_NORMALIZATION_LENGTH = 3.49869654884664


# =========================================================
# DIRECT CALCULATION
# =========================================================

results = run_vfe_vme_calculation(
    potential_file=POTENTIAL_FILE,
    lattice_parameter=LATTICE_PARAMETER,
    composition=COMPOSITION,
    fcc_coordination_file=FCC_COORDINATION_FILE,
    vacancy_environment_file=VACANCY_ENVIRONMENT_FILE,
    transition_state_environment_file=TRANSITION_STATE_ENVIRONMENT_FILE,
    use_alpha=USE_ALPHA,
    alpha_file=ALPHA_FILE,
    environment_normalization_length=ENVIRONMENT_NORMALIZATION_LENGTH,
)


# =========================================================
# PRINT RESULTS
# =========================================================

print_summary(
    results,
    precision=8,
)


# =========================================================
# RESULT TABLES
# =========================================================

fcc_table = results["FCC_statistics"]

vfe_table = results["VFE_site_statistics"]

ts_table = results["TS_site_statistics"]


# =========================================================
# OPTIONAL CSV EXPORTS
# =========================================================

# fcc_table.to_csv(
#     "FCC_statistics.csv",
#     index=False,
# )

# vfe_table.to_csv(
#     "VFE_site_statistics.csv",
#     index=False,
# )

# ts_table.to_csv(
#     "TS_site_statistics.csv",
#     index=False,
# )