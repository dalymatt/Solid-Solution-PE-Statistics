# -*- coding: utf-8 -*-
"""
Reusable workflow for analytical cohesive-energy, VFE/VME, and GPFE
statistics in concentrated FCC alloys.

This module orchestrates the equation implementations in ``Potential.py``,
``Potential_GPFE.py``, and ``defect_energy_functions.py``. Equation references:
cohesive/SRO statistics = vacancy-manuscript Eqs. (1)-(9); VFE/VME = Eqs.
(10)-(14); random-alloy GPFE = CMS-2022 Eqs. (11)-(16).

Public function
---------------
run_energy_calculations()

All calculation switches, fault types, the paths of the defectenvironment files in
``data/coordination``, and the environment normalization length are
defined in this module. Perfect-FCC structure factors are generated from
the lattice parameter and cutoff radius.
"""

from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple, Union
import os

import numpy as np
import pandas as pd


PathLike = Union[str, Path]


# =========================================================
# DEFAULT CALCULATION SETTINGS
# =========================================================

CALCULATE_COHESIVE = True
CALCULATE_VFE_VME = True
CALCULATE_GPFE = True

FAULT_TYPES = (
    "USF",
    "ISF",
    "UTF1",
    "ESF",
    "UTF2",
    "TF",
)

ENVIRONMENT_NORMALIZATION_LENGTH = 3.49869654884664


# =========================================================
# INPUT VALIDATION
# =========================================================

def _validate_inputs(
    *,
    potential_file: PathLike,
    composition: Sequence[float],
    mode: str,
    alpha_file: Optional[PathLike],
) -> Tuple[np.ndarray, bool, Optional[Path]]:
    """Validate common material inputs and select Random or SRO mode."""

    potential_path = Path(potential_file).resolve()
    composition_array = np.asarray(composition, dtype=float)

    if composition_array.ndim != 1:
        raise ValueError("COMPOSITION must be a one-dimensional array.")

    if composition_array.size == 0:
        raise ValueError("COMPOSITION cannot be empty.")

    if not np.all(np.isfinite(composition_array)):
        raise ValueError("COMPOSITION contains NaN or infinite values.")

    if np.any(composition_array < 0.0):
        raise ValueError("COMPOSITION cannot contain negative values.")

    if not np.isclose(
        composition_array.sum(),
        1.0,
        atol=1.0e-6,
    ):
        raise ValueError(
            "COMPOSITION must sum to 1. "
            f"Current sum = {composition_array.sum():.10f}"
        )

    if not potential_path.is_file():
        raise FileNotFoundError(
            f"EAM potential file was not found:\n{potential_path}"
        )

    normalized_mode = mode.strip().lower()

    if normalized_mode == "random":
        return composition_array, False, None

    if normalized_mode == "sro":
        if alpha_file is None:
            raise ValueError(
                "ALPHA_FILE must be supplied when MODE = 'SRO'."
            )

        alpha_path = Path(alpha_file).resolve()

        if not alpha_path.is_file():
            raise FileNotFoundError(
                "Warren-Cowley parameter file was not found:\n"
                f"{alpha_path}"
            )

        return composition_array, True, alpha_path

    raise ValueError('MODE must be either "Random" or "SRO".')

# =========================================================
# FCC COORDINATION AND SRO PARAMETERS
# =========================================================

def _generate_fcc_coordination(
    *,
    lattice_parameter: float,
    cutoff_radius: float,
):
    """Generate perfect-FCC RDF and coordination relations."""

    import rdf_coord as rc

    return rc.rdf_coord_fcc(
        lattice_parameter,
        cutoff_radius,
    )


def _load_alpha(
    *,
    use_alpha: bool,
    alpha_file: Optional[PathLike],
    fcc_coordination: np.ndarray,
    composition: np.ndarray,
) -> np.ndarray:
    """Load Warren-Cowley parameters or create zeros for a random alloy."""

    expected_shape = (
        fcc_coordination.shape[0],
        composition.size,
        composition.size,
    )

    if not use_alpha:
        return np.zeros(expected_shape, dtype=float)

    alpha = np.asarray(
        np.load(alpha_file),
        dtype=float,
    )

    if alpha.shape != expected_shape:
        raise ValueError(
            "Warren-Cowley array shape mismatch.\n"
            f"Expected: {expected_shape}\n"
            f"Received: {alpha.shape}\n"
            "Check CUTOFF_RADIUS, COMPOSITION, and ALPHA_FILE."
        )

    if not np.all(np.isfinite(alpha)):
        raise ValueError(
            "The Warren-Cowley array contains NaN or infinite values."
        )

    return alpha


# =========================================================
# COHESIVE-ENERGY STATISTICS
# =========================================================

def _calculate_cohesive_statistics(
    *,
    potential_file: PathLike,
    composition: np.ndarray,
    fcc_coordination: np.ndarray,
    alpha: np.ndarray,
) -> Dict[str, Any]:
    """Calculate FCC cohesive-energy statistics.

    Delegates to ``Potential.potential_stats`` (vacancy manuscript
    Eqs. (1)-(9); random-alloy limit: CMS-2022 Eqs. (3)-(10)).
    """

    import Potential as pot

    rrange, rhorange, rho, Fr, Pp = pot.potential_read(
        str(potential_file)
    )

    table, covariance, element_energies = pot.potential_stats(
        rrange,
        rhorange,
        rho,
        Fr,
        Pp,
        composition,
        fcc_coordination,
        alpha,
    )
    
    if not np.all(np.isfinite(table["E"].to_numpy(dtype=float))):
        raise ValueError("Non-finite cohesive-energy statistics encountered; check variance/covariance inputs.")

    return {
        "table": table,
        "covariance": covariance,
        "element_energies": element_energies,
    }


def _print_cohesive_summary(
    *,
    table: pd.DataFrame,
    mode: str,
) -> None:
    """Print FCC cohesive-energy statistics."""

    print()
    print("=" * 66)
    print(
        f"FCC COHESIVE-ENERGY STATISTICS: "
        f"{mode.upper()} SYSTEM"
    )
    print("=" * 66)
    print(table)
    print("=" * 66)


# =========================================================
# VFE/VME STATISTICS
# =========================================================

def _calculate_vfe_vme_statistics(
    *,
    potential_file: PathLike,
    lattice_parameter: float,
    composition: np.ndarray,
    fcc_coordination: np.ndarray,
    vacancy_environment_file: PathLike,
    transition_state_environment_file: PathLike,
    use_alpha: bool,
    alpha_file: Optional[PathLike],
) -> Dict[str, Any]:
    """Calculate VFE/VME statistics from pre-exported environments.

    Mean energies follow vacancy-manuscript Eqs. (10)-(11); standard
    deviations follow Eqs. (12)-(14).
    """

    from defect_energy_functions import run_vfe_vme_calculation

    required_files = (        
        Path(vacancy_environment_file),
        Path(transition_state_environment_file),
    )

    for required_file in required_files:
        if not required_file.is_file():
            raise FileNotFoundError(
                "Required VFE/VME coordination file was not found:\n"
                f"{required_file}"
            )

    return run_vfe_vme_calculation(
        potential_file=potential_file,
        lattice_parameter=lattice_parameter,
        composition=composition,
        fcc_coordination=fcc_coordination,
        vacancy_environment_file=vacancy_environment_file,
        transition_state_environment_file=(
            transition_state_environment_file
        ),
        use_alpha=use_alpha,
        alpha_file=alpha_file,
        environment_normalization_length=(
            ENVIRONMENT_NORMALIZATION_LENGTH
        ),
    )


def _print_defect_summary(
    results: Dict[str, Any],
    precision: int = 8,
) -> None:
    """Print VFE, transition-state, and VME statistics."""

    print()
    print("=" * 60)
    print("DEFECT ENERGY SUMMARY")
    print("=" * 60)

    print(
        f"VFE average                 = "
        f"{results['E_VFE']:.{precision}f} eV"
    )
    print(
        f"VFE stdev                   = "
        f"{results['Sig_VFE']:.{precision}f} eV"
    )

    print("-" * 60)

    print(
        f"TS excess energy average    = "
        f"{results['E_TS']:.{precision}f} eV"
    )
    print(
        f"TS excess energy stdev      = "
        f"{results['Sig_TS']:.{precision}f} eV"
    )

    print("-" * 60)

    print(
        f"VME average                 = "
        f"{results['E_VME']:.{precision}f} eV"
    )
    print(
        f"VME stdev                   = "
        f"{results['Sig_VME']:.{precision}f} eV"
    )

    print("=" * 60)


# =========================================================
# GPFE STATISTICS
# =========================================================

def _calculate_gpfe_statistics(
    *,
    root: PathLike,
    potential_file: PathLike,
    lattice_parameter: float,
    cutoff_radius: float,
    composition: np.ndarray,
    fcc_coordination: np.ndarray,
) -> Dict[str, Any]:
    """Calculate random-alloy GPFE statistics.

    Fault means and fluctuations follow CMS-2022 Eqs. (11)-(16).
    """

    import Potential_GPFE as pot_gpfe
    import rdf_coord as rc

    root_path = Path(root).resolve()
    coordination_directory = (
        root_path / "data" / "coordination"
    )

    if not coordination_directory.is_dir():
        raise FileNotFoundError(
            "The coordination-data directory was not found:\n"
            f"{coordination_directory}"
        )

    rrange, rhorange, rho, Fr, Pp = (
        pot_gpfe.potential_read(
            str(potential_file)
        )
    )

    (
        fcc_table,
        fcc_covariance,
        fcc_element_energies,
    ) = pot_gpfe.potential_stats(
        rrange,
        rhorange,
        rho,
        Fr,
        Pp,
        composition,
        fcc_coordination,
    )

    fault_results: Dict[str, pd.DataFrame] = {}
    fault_covariances: Dict[str, float] = {}
    fault_coordination: Dict[str, Any] = {}

    for fault_type in FAULT_TYPES:

        fault_file = (
            coordination_directory
            / f"cn_{fault_type}.pkl"
        )

        if not fault_file.is_file():
            raise FileNotFoundError(
                "Fault coordination file was not found:\n"
                f"{fault_file}"
            )

        _, cn_fault = rc.rdf_coord_fault(
            lattice_parameter,
            cutoff_radius,
            fcc_coordination,
            fault_type,
        )

        table, covariance = pot_gpfe.potential_stats_fault(
            rrange,
            rhorange,
            rho,
            Fr,
            Pp,
            composition,
            fcc_coordination,
            cn_fault,
            fcc_table,
            fcc_element_energies,
            fault_type,
        )

        fault_results[fault_type] = table
        fault_covariances[fault_type] = float(covariance)
        fault_coordination[fault_type] = cn_fault

    summary_table = pd.concat(
        [
            fault_results[fault_type]
            for fault_type in FAULT_TYPES
        ],
        axis=1,
    )
    
    # Convert the internal GPFE quantity to the conventional surface-energy
    # unit.  For FCC {111}, rho_111 = 4/(sqrt(3)*a^2) atoms/A^2 and
    # 1 eV/A^2 = 16021.76634 mJ/m^2.
    rho_111 = 4.0 / (np.sqrt(3.0) * lattice_parameter**2)
    gpfe_conversion = rho_111 * 16021.76634
    summary_mj_m2 = summary_table * gpfe_conversion
    
    

    return {
        "summary": summary_table,
        "summary_mJ_m2": summary_mj_m2,
        "conversion_eVatom_to_mJm2": gpfe_conversion,
        "results": fault_results,
        "covariances": fault_covariances,
        "coordination": fault_coordination,
        "fcc_table": fcc_table,
        "fcc_covariance": fcc_covariance,
        "fcc_element_energies": fcc_element_energies,
    }


def _print_gpfe_summary(
    table: pd.DataFrame,
    table_mj_m2: Optional[pd.DataFrame] = None,
) -> None:
    """Print the combined GPFE table."""

    print()
    print("=" * 66)
    print("GENERALIZED PLANAR FAULT ENERGY STATISTICS")
    print("=" * 66)
    print(table)
    if table_mj_m2 is not None:
        print("-" * 66)
        print("GENERALIZED PLANAR FAULT ENERGY STATISTICS (mJ/m^2)")
        print("-" * 66)
        print(table_mj_m2)
    print("=" * 66)


# =========================================================
# PUBLIC WORKFLOW
# =========================================================

def run_energy_calculations(
    *,
    root: PathLike,
    potential_file: PathLike,
    lattice_parameter: float,
    cutoff_radius: float,
    composition: Sequence[float],
    mode: str,
    alpha_file: Optional[PathLike] = None,
    precision: int = 8,
) -> Dict[str, Any]:
    """
    Run all enabled energy calculations.

    Standard paths are determined automatically from ROOT:
    
    data/coordination/cn_vac.pkl
    data/coordination/cn_TS.pkl

    Returns
    -------
    dict
        Contains the complete cohesive, defect, and GPFE results, along
        with convenient table aliases.
    """

    root_path = Path(root).resolve()

    (
        composition_array,
        use_alpha,
        active_alpha_file,
    ) = _validate_inputs(
        potential_file=potential_file,
        composition=composition,
        mode=mode,
        alpha_file=alpha_file,
    )

    _, fcc_coordination = _generate_fcc_coordination(
        lattice_parameter=lattice_parameter,
        cutoff_radius=cutoff_radius,
    )
    # Dimensionless copy for the defect path; run_vfe_vme_calculation
    # multiplies the first column by the lattice parameter.
    fcc_coordination_normalized = fcc_coordination.copy()
    fcc_coordination_normalized[:, 0] /= lattice_parameter

    alpha = _load_alpha(
        use_alpha=use_alpha,
        alpha_file=active_alpha_file,
        fcc_coordination=fcc_coordination,
        composition=composition_array,
    )

    vacancy_environment_file = (
        root_path / "data" / "coordination" / "cn_vac.pkl"
    )

    transition_state_environment_file = (
        root_path / "data" / "coordination" / "cn_TS.pkl"
    )

    output: Dict[str, Any] = {
        "mode": mode,
        "use_alpha": use_alpha,
        "alpha": alpha,
        "fcc_coordination": fcc_coordination,
        "cohesive": None,
        "defect": None,
        "gpfe": None,
        "cohesive_table": None,
        "vfe_table": None,
        "ts_table": None,
        "gpfe_table": None,
    }

    if CALCULATE_COHESIVE:

        cohesive = _calculate_cohesive_statistics(
            potential_file=potential_file,
            composition=composition_array,
            fcc_coordination=fcc_coordination,
            alpha=alpha,
        )

        _print_cohesive_summary(
            table=cohesive["table"],
            mode=mode,
        )

        output["cohesive"] = cohesive
        output["cohesive_table"] = cohesive["table"]

    if CALCULATE_VFE_VME:

        defect = _calculate_vfe_vme_statistics(
            potential_file=potential_file,
            lattice_parameter=lattice_parameter,
            composition=composition_array,
            fcc_coordination=fcc_coordination_normalized,
            vacancy_environment_file=(
                vacancy_environment_file
            ),
            transition_state_environment_file=(
                transition_state_environment_file
            ),
            use_alpha=use_alpha,
            alpha_file=active_alpha_file,
        )

        _print_defect_summary(
            defect,
            precision=precision,
        )

        output["defect"] = defect
        output["vfe_table"] = defect[
            "VFE_site_statistics"
        ]
        output["ts_table"] = defect[
            "TS_site_statistics"
        ]

    if CALCULATE_GPFE and not use_alpha:

        gpfe = _calculate_gpfe_statistics(
            root=root_path,
            potential_file=potential_file,
            lattice_parameter=lattice_parameter,
            cutoff_radius=cutoff_radius,
            composition=composition_array,
            fcc_coordination=fcc_coordination,
        )

        _print_gpfe_summary(
            gpfe["summary"],
            gpfe["summary_mJ_m2"],
        )

        output["gpfe"] = gpfe
        output["gpfe_table"] = gpfe["summary"]
        output["gpfe_table_mJ_m2"] = gpfe["summary_mJ_m2"]

    elif CALCULATE_GPFE and use_alpha:

        print()
        print("=" * 66)
        print("GPFE CALCULATION SKIPPED")
        print("=" * 66)
        print(
            "MODE = 'SRO'. GPFE calculations are currently "
            "implemented for random alloys only."
        )
        print("=" * 66)

    return output