# -*- coding: utf-8 -*-
"""
Example input script for the analytical calculation of cohesive-energy,
vacancy-formation-energy (VFE), vacancy-migration-energy (VME), and
generalized planar-fault-energy (GPFE) statistics in concentrated
face-centered-cubic (FCC) solid solutions.

This file is intended to be executed directly. It contains only the
material-specific inputs. The reusable analytical equations and numerical
procedures are implemented in:

    src/energy_workflow.py
    src/defect_energy_functions.py
    src/Potential.py
    src/Potential_GPFE.py

Analytical basis
----------------
1. R. Jagatramka, C. Wang, and M. Daly,
   "An analytical method to quantify the statistics of energy landscapes
   in random solid solutions,"
   Computational Materials Science 214 (2022) 111763.
   https://doi.org/10.1016/j.commatsci.2022.111763

2. A. Baski, R. Jagatramka, and M. Daly,
   "A mechanistic model for vacancy 
       energetics in concentrated solid solutions with short-range order"

The workflow evaluates statistical energy landscapes using an EAM/alloy
potential, FCC coordination structure factors, alloy composition, and,
when requested, Warren-Cowley short-range-order parameters.

"""

from pathlib import Path
import sys

import numpy as np


# =============================================================================
# REPOSITORY PATHS
# =============================================================================

# ROOT is the directory containing this example script.
#
# Expected repository organization:
#
# repository_root/
# |-- example_defect_energy.py
# |-- src/
# |   |-- energy_workflow.py
# |   |-- defect_energy_functions.py
# |   |-- Potential.py
# |   `-- Potential_GPFE.py
# |-- potentials/
# |-- data/
# `-- ...
ROOT = Path(__file__).resolve().parent

# Directory containing the reusable calculation modules.
SRC_DIR = ROOT / "src"

# Add src/ to Python's module-search path so that energy_workflow.py can be
# imported when this script is executed directly from the repository root.
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# Main public workflow function. It coordinates the cohesive-energy,
# VFE/VME, and GPFE calculations using the inputs defined below.
from energy_workflow import run_energy_calculations


# =============================================================================
# MATERIAL-SPECIFIC INPUTS
# =============================================================================

# -----------------------------------------------------------------------------
# EAM/ALLOY POTENTIAL FILE
# -----------------------------------------------------------------------------
#
# Select the setfl-style EAM/alloy potential that describes the alloy system.
#
# IMPORTANT:
# The element order stored in the potential file determines the required order
# of the entries in COMPOSITION.
#
# Active example:
#     Fe-Ni-Cr EAM/alloy potential
POTENTIAL_FILE = ROOT / "potentials" / "FeNiCr.eam.alloy"

# Alternative example for an equimolar Fe-Ni-Cr-Co-Cu alloy:
#POTENTIAL_FILE = ROOT / "potentials" / "FeNiCrCoCu-with-ZBL.eam.alloy" 


# -----------------------------------------------------------------------------
# FCC LATTICE PARAMETER
# -----------------------------------------------------------------------------
#
# Equilibrium or prescribed FCC lattice parameter in angstroms.
#
# The normalized coordination structure factors are multiplied by this value
# to recover physical interatomic distances before the EAM functions are
# interpolated.
#
# Use a lattice parameter consistent with:
#     1. the selected alloy composition,
#     2. the selected random/SRO state, and
#     3. the selected EAM potential.
#LATTICE_PARAMETER = 3.5225        # Fe0.33Ni0.33Cr0.34, random alloy

# Additional examples:
LATTICE_PARAMETER = 3.51036     # Fe0.73Ni0.08Cr0.19, random alloy


#LATTICE_PARAMETER = 3.53073      # Equimolar NiCrCo
#LATTICE_PARAMETER = 3.546        # Equimolar NiCrCoCu
#LATTICE_PARAMETER = 3.54939      # Equimolar FeNiCrCoCu


# -----------------------------------------------------------------------------
# NEIGHBOR CUTOFF RADIUS
# -----------------------------------------------------------------------------
#
# Radial cutoff in angstroms used to construct the FCC radial-distribution and
# coordination-shell relations.
#
# The cutoff must be consistent with the selected potential and with the
# coordination data used by the analytical defect-energy calculations.
CUTOFF_RADIUS = 5.6               # FeNiCr potential

# Alternative example:
#CUTOFF_RADIUS = 5.80375         # FeNiCrCoCu potential


# -----------------------------------------------------------------------------
# ALLOY COMPOSITION
# -----------------------------------------------------------------------------
#
# Atomic fractions of the alloy components.
#
# Requirements:
#     1. COMPOSITION must be one-dimensional.
#     2. Every value must be non-negative.
#     3. The entries must sum to 1.0.
#     4. The element order must match the order in POTENTIAL_FILE.
#
# For the FeNiCr potential, the expected order is:
#
#     [Fe, Ni, Cr]
#
# Thus, the active input corresponds to:
#
#     Fe0.33 Ni0.33 Cr0.34
#COMPOSITION = np.array([0.33, 0.33, 0.34], dtype=float)

# Additional examples:
COMPOSITION = np.array([0.73, 0.08, 0.19], dtype=float)
#
# For a five-component potential whose element order is
# [Fe, Ni, Cr, Co, Cu]:
#COMPOSITION = np.array([0.2, 0.2, 0.2, 0.2, 0.2], dtype=float)
#COMPOSITION = np.array([0.0, 0.333, 0.333, 0.334, 0.0], dtype=float)
#COMPOSITION = np.array([0.0, 0.25, 0.25, 0.25, 0.25], dtype=float)


# -----------------------------------------------------------------------------
# CHEMICAL-ARRANGEMENT MODE
# -----------------------------------------------------------------------------
#
# Select exactly one of the following:
#
#     MODE = "Random"
#         All Warren-Cowley parameters are set to zero:
#
#             alpha_ij^(zeta) = 0
#
#         This represents statistically random occupation of the lattice.
#
#     MODE = "SRO"
#         Warren-Cowley parameters are loaded from ALPHA_FILE and used to
#         modify the conditional neighbor probabilities in each coordination
#         shell.
#
# The Warren-Cowley parameter follows the convention
#
#     alpha_ij^(zeta) = 1 - P_ij^(zeta) / c_j,
#
# where P_ij^(zeta) is the conditional probability of finding species j
# around species i in shell zeta, and c_j is the global concentration of j.
MODE = "Random"
#MODE = "SRO"


# -----------------------------------------------------------------------------
# WARREN-COWLEY SRO PARAMETER FILE
# -----------------------------------------------------------------------------
#
# NumPy .npy file containing the Warren-Cowley parameter array.
#
# Expected array shape:
#
#     (number_of_coordination_shells,
#      number_of_elements,
#      number_of_elements)
#
# This file is used only when MODE = "SRO". It may remain defined while
# MODE = "Random"; in that case, the workflow ignores it and uses alpha = 0.


ALPHA_FILE = (ROOT / "data" / "sro" / "304SS" / "alpha_Fe73Ni8Cr19_alpha_p0p10.npy")
#ALPHA_FILE = (ROOT / "data" / "sro" / "Quinary"  / "alpha_Fe20Ni20Cr20Co20Cu20_alpha_P0p15.npy")


# =============================================================================
# RUN ANALYTICAL CALCULATIONS
# =============================================================================
#
# run_energy_calculations() coordinates the following operations:
#
# 1. Validate the potential path, composition, and random/SRO selection.
#
# 2. Generate the perfect-FCC coordination relations.
#
# 3. Read the EAM/alloy functions:
#       - elemental electron-density functions,
#       - embedding-energy functions, and
#       - pair-interaction functions.
#
# 4. Calculate FCC cohesive-energy statistics using the reparameterized EAM
#    framework. The underlying energy decomposition corresponds to the
#    embedding and pair contributions described in the analytical papers.
#
# 5. Calculate vacancy-formation-energy statistics from the excess energy of
#    vacancy-containing coordination environments relative to perfect FCC.
#
# 6. Calculate transition-state excess-energy statistics and obtain the
#    vacancy migration-energy distribution from the difference between the
#    transition-state and vacancy-state energy landscapes.
#
# 7. Calculate GPFE statistics for the enabled planar-fault configurations.
#
# The returned object is a dictionary containing both detailed tables and
# scalar summary quantities.
results = run_energy_calculations(
    root=ROOT,
    potential_file=POTENTIAL_FILE,
    lattice_parameter=LATTICE_PARAMETER,
    cutoff_radius=CUTOFF_RADIUS,
    composition=COMPOSITION,
    mode=MODE,
    alpha_file=ALPHA_FILE,
)


# =============================================================================
# OPTIONAL RESULT VARIABLES
# =============================================================================
#
# The workflow already prints the enabled calculation summaries. The aliases
# below are provided for interactive analysis, plotting, exporting, or further
# post-processing after this script has run.
#
# Each variable remains a reference to the corresponding object stored in the
# results dictionary.


# FCC cohesive-energy statistics.
#
# Typical rows:
#     Mean
#     Std
#
# Typical columns include:
#     rho  : local electron-density contribution
#     F    : embedding-energy contribution
#     Pp   : pair-interaction contribution
#     E    : total per-atom energy
cohesive_table = results["cohesive_table"]


# Site-resolved vacancy-formation-energy information.
#
# This table contains statistics associated with each unique vacancy
# coordination environment used to assemble the total VFE distribution.
vfe_table = results["vfe_table"]


# Site-resolved transition-state excess-energy information.
#
# This table contains the statistics of the unique saddle-point coordination
# environments used in the analytical VME calculation.
ts_table = results["ts_table"]


# Generalized planar-fault-energy statistics.
#
# Depending on the calculation settings in energy_workflow.py, this table may
# contain results for unstable stacking fault, intrinsic stacking fault,
# unstable twinning fault, extrinsic stacking fault, and related fault states.
gpfe_table = results["gpfe_table"]


# =============================================================================
# ADDITIONAL SCALAR OUTPUTS
# =============================================================================
#
# The results dictionary may also contain scalar quantities such as:
#
#     results["E_VFE"]    : average vacancy formation energy
#     results["Sig_VFE"]  : VFE standard deviation
#     results["E_TS"]     : average transition-state excess energy
#     results["Sig_TS"]   : transition-state standard deviation
#     results["E_VME"]    : average vacancy migration energy
#     results["Sig_VME"]  : VME standard deviation
#
# These entries can be accessed directly when needed, for example:
#
# average_vfe = results["E_VFE"]
# std_vfe = results["Sig_VFE"]
# average_vme = results["E_VME"]
# std_vme = results["Sig_VME"]
