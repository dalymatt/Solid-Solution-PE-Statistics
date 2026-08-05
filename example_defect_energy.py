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

from pathlib import Path
import sys

import numpy as np

# Repository paths: users can run this file directly from any working directory.
ROOT = Path(__file__).resolve().parent
SRC_DIR = ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

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
POTENTIAL_FILE = ROOT / "potentials" / "FeNiCr.eam.alloy"
#POTENTIAL_FILE = ROOT / "potentials" / "FeNiCrCoCu-with-ZBL.eam.alloy"

# Lattice parameter in angstrom
LATTICE_PARAMETER = 3.51036        # Fe73Ni8Cr19 random
# LATTICE_PARAMETER = 3.54939      # Equimolar FeNiCrCoCu
# LATTICE_PARAMETER = 3.4986       # Pure Fe in Fe-Ni-Cr potential
# LATTICE_PARAMETER = 3.53073      # Equimolar NiCrCo
# LATTICE_PARAMETER = 3.546        # Equimolar NiCrCoCu
# LATTICE_PARAMETER = 3.50931      # Fe73Ni8Cr19, alpha = +0.05
# LATTICE_PARAMETER = 3.5150       # Fe73Ni8Cr19, alpha = -0.05

# Cutoff radius in angstrom
CUTOFF_RADIUS = 5.6                         # Fe-Ni-Cr
# CUTOFF_RADIUS = 5.80375                   # Fe-Ni-Cr-Co-Cu

# Composition must follow the element order in the potential file.
COMPOSITION = np.array([0.73, 0.08, 0.19])       # FeNiCr
# COMPOSITION = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
# COMPOSITION = np.array([0.0, 0.333, 0.333, 0.334, 0.0])
# COMPOSITION = np.array([0.0, 0.25, 0.25, 0.25, 0.25])
# COMPOSITION = np.array([1.0, 0.0, 0.0])


# =========================================================
# PRE-EXPORTED COORDINATION FILES
# =========================================================

# Normalized perfect-FCC coordination structure
FCC_COORDINATION_FILE = ROOT / "data" / "coordination" / "cn_FCC.pkl"

# Vacancy environments
VACANCY_ENVIRONMENT_FILE = ROOT / "data" / "coordination" / "cn_vac.pkl"

# Transition-state environments
TRANSITION_STATE_ENVIRONMENT_FILE = ROOT / "data" / "coordination" / "cn_TS.pkl."


# =========================================================
# ORDERING CONDITION
# =========================================================

# Select exactly one mode: "Random" or "SRO".
MODE = "SRO"

if MODE == "Random":
    USE_ALPHA = False
    ALPHA_FILE = None
elif MODE == "SRO":
    USE_ALPHA = True
    ALPHA_FILE = (
        ROOT / "data" / "sro"
        / "alpha_FeNiCr_SS_minus_point05.npy"
    )
else:
    raise ValueError('MODE must be either "Random" or "SRO".')


# =========================================================
# DEFECT-ENVIRONMENT NORMALIZATION
# =========================================================

ENVIRONMENT_NORMALIZATION_LENGTH = 3.49869654884664 # Normalized from Pure Fe Lattice Parameter


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