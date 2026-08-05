# -*- coding: utf-8 -*-
"""
Direct-execution example for:

1. FCC cohesive energy
2. Vacancy formation energy (VFE)
3. Vacancy migration energy (VME)
4. Generalized planar fault energies (GPFEs)

Current implementation
----------------------
Random system:
    FCC + VFE + VME + GPFE

SRO system:
    FCC + VFE + VME

GPFE is automatically skipped when USE_ALPHA = True.
"""

import numpy as np

from defect_gpfe_functions import (
    calculate_random_gpfe,
    print_gpfe_summary,
    print_summary,
    run_vfe_vme_calculation,
)


# =========================================================
# SYSTEM INPUTS
# =========================================================

# EAM/alloy potential file
POTENTIAL_FILE = "FeNiCr.eam.alloy"

# Lattice parameter in angstrom
LATTICE_PARAMETER = 3.51036

# Cutoff radius in angstrom
CUTOFF_RADIUS = 5.6

# Composition must follow the element order in the potential
COMPOSITION = np.array([
    0.73,
    0.08,
    0.19,
])


# =========================================================
# PRE-EXPORTED COORDINATION FILES
# =========================================================

# Normalized perfect-FCC coordination structure
FCC_COORDINATION_FILE = "cn_FCC.pkl"

# Vacancy environments
VACANCY_ENVIRONMENT_FILE = "cn_vac.pkl"

# Transition-state environments
TRANSITION_STATE_ENVIRONMENT_FILE = "TS_rel.pkl"


# =========================================================
# ORDERING CONDITION
# =========================================================

# ---------------------------------------------------------
# Random system
# ---------------------------------------------------------

#USE_ALPHA = False
#ALPHA_FILE = None


# ---------------------------------------------------------
# SRO system example
# ---------------------------------------------------------

USE_ALPHA = True
ALPHA_FILE = "alpha_FeNiCr_SS_minus_point05.npy"


# =========================================================
# DEFECT-ENVIRONMENT NORMALIZATION
# =========================================================

ENVIRONMENT_NORMALIZATION_LENGTH = 3.49869654884664


# =========================================================
# SELECT CALCULATIONS
# =========================================================

# VFE and VME are calculated for both random and SRO systems
CALCULATE_VFE_VME = True

# GPFE is automatically calculated only for random systems
CALCULATE_GPFE = not USE_ALPHA


# =========================================================
# DIRECT VFE/VME CALCULATION
# =========================================================

if CALCULATE_VFE_VME:

    defect_results = run_vfe_vme_calculation(
        potential_file=POTENTIAL_FILE,
        lattice_parameter=LATTICE_PARAMETER,
        composition=COMPOSITION,
        fcc_coordination_file=FCC_COORDINATION_FILE,
        vacancy_environment_file=VACANCY_ENVIRONMENT_FILE,
        transition_state_environment_file=(
            TRANSITION_STATE_ENVIRONMENT_FILE
        ),
        use_alpha=USE_ALPHA,
        alpha_file=ALPHA_FILE,
        environment_normalization_length=(
            ENVIRONMENT_NORMALIZATION_LENGTH
        ),
    )

    print_summary(
        defect_results,
        precision=8,
    )

    # Store calculated tables
    fcc_table = defect_results["FCC_statistics"]

    vfe_table = defect_results[
        "VFE_site_statistics"
    ]

    ts_table = defect_results[
        "TS_site_statistics"
    ]

    # Store individual values
    E_FCC = defect_results["E_FCC"]
    Sig_FCC = defect_results["Sig_FCC"]

    E_VFE = defect_results["E_VFE"]
    Sig_VFE = defect_results["Sig_VFE"]

    E_TS = defect_results["E_TS"]
    Sig_TS = defect_results["Sig_TS"]

    E_VME = defect_results["E_VME"]
    Sig_VME = defect_results["Sig_VME"]


# =========================================================
# DIRECT RANDOM-ALLOY GPFE CALCULATION
# =========================================================

if CALCULATE_GPFE:

    gpfe_results = calculate_random_gpfe(
        potential_file=POTENTIAL_FILE,
        lattice_parameter=LATTICE_PARAMETER,
        cutoff_radius=CUTOFF_RADIUS,
        composition=COMPOSITION,
        fcc_coordination_file=(
            FCC_COORDINATION_FILE
        ),
        fault_types=(
            "USF",
            "ISF",
            "UTF1",
            "ESF",
            "UTF2",
            "TF",
        ),
    )

    print_gpfe_summary(
        gpfe_results,
        precision=8,
    )

    # Combined GPFE summary table
    gpfe_table = gpfe_results[
        "GPFE_summary"
    ]

    # Individual fault-energy tables
    form_E_USF = gpfe_results[
        "GPFE_results"
    ]["USF"]

    form_E_ISF = gpfe_results[
        "GPFE_results"
    ]["ISF"]

    form_E_UTF1 = gpfe_results[
        "GPFE_results"
    ]["UTF1"]

    form_E_ESF = gpfe_results[
        "GPFE_results"
    ]["ESF"]

    form_E_UTF2 = gpfe_results[
        "GPFE_results"
    ]["UTF2"]

    form_E_TF = gpfe_results[
        "GPFE_results"
    ]["TF"]


# =========================================================
# AUTOMATIC GPFE SKIP FOR SRO
# =========================================================

else:

    gpfe_results = None
    gpfe_table = None

    form_E_USF = None
    form_E_ISF = None
    form_E_UTF1 = None
    form_E_ESF = None
    form_E_UTF2 = None
    form_E_TF = None

    print()
    print("=" * 58)
    print("GPFE CALCULATION SKIPPED")
    print("=" * 58)
    print(
        "USE_ALPHA = True, so this is an SRO calculation."
    )
    print(
        "Only FCC, VFE, and VME are currently calculated "
        "for SRO systems."
    )
    print(
        "GPFE for SRO systems will be included in a future study."
    )
    print("=" * 58)


# =========================================================
# OPTIONAL CSV EXPORTS
# =========================================================

# ---------------------------------------------------------
# FCC, VFE, and TS tables
# ---------------------------------------------------------

# if CALCULATE_VFE_VME:
#
#     fcc_table.to_csv(
#         "FCC_statistics.csv"
#     )
#
#     vfe_table.to_csv(
#         "VFE_site_statistics.csv",
#         index=False,
#     )
#
#     ts_table.to_csv(
#         "TS_site_statistics.csv",
#         index=False,
#     )


# ---------------------------------------------------------
# GPFE summary table
# ---------------------------------------------------------

# if CALCULATE_GPFE:
#
#     gpfe_table.to_csv(
#         "GPFE_summary.csv",
#         index=False,
#     )