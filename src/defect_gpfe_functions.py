# -*- coding: utf-8 -*-
"""
Reusable functions for analytical FCC, VFE, and VME calculations.

The perfect-FCC coordination/structure-factor array is loaded directly
from a pre-exported pickle file.
"""

from __future__ import annotations

import copy
import pickle
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple, Union

import numpy as np
import pandas as pd

import Potential as pot_defect
import Potential_GPFE as pot_gpfe
import rdf_coord as rc


PathLike = Union[str, Path]


def load_pickle(filename: PathLike) -> Any:
    """Load and return an object stored in a pickle file."""
    with open(filename, "rb") as file_handle:
        return pickle.load(file_handle)


def load_fcc_coordination(fcc_coordination_file: PathLike) -> np.ndarray:
    """
    Load the normalized perfect-FCC coordination/structure-factor array.

    The first column must already be dimensionless. It will later be
    multiplied by the requested lattice parameter.

    Parameters
    ----------
    fcc_coordination_file
        Pickle file containing ``cn_FCC``.

    Returns
    -------
    numpy.ndarray
        FCC coordination array.
    """
    cn_fcc = np.asarray(
        load_pickle(fcc_coordination_file),
        dtype=float,
    )

    if cn_fcc.ndim != 2:
        raise ValueError(
            "The FCC coordination array must be two-dimensional. "
            f"Received shape {cn_fcc.shape}."
        )

    if cn_fcc.shape[1] < 2:
        raise ValueError(
            "The FCC coordination array must contain at least two columns."
        )

    if not np.all(np.isfinite(cn_fcc)):
        raise ValueError(
            "The FCC coordination array contains NaN or infinite values."
        )

    return cn_fcc


def load_alpha_parameters(
    alpha_file: PathLike | None,
    use_alpha: bool,
    number_of_shells: int,
    number_of_elements: int,
) -> np.ndarray:
    """
    Load Warren-Cowley parameters or create zeros for a random alloy.
    """
    expected_shape = (
        number_of_shells,
        number_of_elements,
        number_of_elements,
    )

    if not use_alpha:
        return np.zeros(expected_shape, dtype=float)

    if alpha_file is None:
        raise ValueError(
            "alpha_file must be supplied when use_alpha=True."
        )

    alpha = np.asarray(np.load(alpha_file), dtype=float)

    if alpha.shape != expected_shape:
        raise ValueError(
            "Alpha shape mismatch: "
            f"expected {expected_shape}, got {alpha.shape}."
        )

    return alpha


def load_coordination_environments(
    pickle_file: PathLike,
    normalization_length: float,
) -> List[np.ndarray]:
    """
    Load and normalize defect coordination environments.

    The first column of each environment is divided by
    ``normalization_length``.
    """
    raw_environments = load_pickle(pickle_file)

    normalized_environments: List[np.ndarray] = []

    for environment in raw_environments:
        normalized = np.asarray(
            copy.deepcopy(environment),
            dtype=float,
        )

        if normalized.ndim != 2:
            raise ValueError(
                f"An environment in {pickle_file} is not two-dimensional."
            )

        normalized[:, 0] /= normalization_length
        normalized_environments.append(normalized)

    return normalized_environments


def get_unique_environments(
    environments: Sequence[np.ndarray],
) -> Tuple[List[np.ndarray], np.ndarray]:
    """Return unique environments and their occurrence frequencies."""
    counts = Counter(
        tuple(map(tuple, np.asarray(environment)))
        for environment in environments
    )

    unique_environments = [
        np.asarray(key, dtype=float)
        for key in counts.keys()
    ]

    frequencies = np.asarray(
        list(counts.values()),
        dtype=int,
    )

    return unique_environments, frequencies


def calculate_fcc_reference(
    potential_file: PathLike,
    composition: np.ndarray,
    lattice_parameter: float,
    cn_fcc: np.ndarray,
    alpha: np.ndarray,
) -> Dict[str, Any]:
    """
    Read the EAM/alloy potential and calculate perfect-FCC statistics.
    """
    rrange, rhorange, rho, Fr, Pp = pot_defect.potential_read(
        str(potential_file)
    )

    cn_perfect = copy.deepcopy(cn_fcc)
    cn_perfect[:, 0] *= lattice_parameter

    coh_fcc, _, element_energy_fcc = pot_defect.potential_stats(
        rrange,
        rhorange,
        rho,
        Fr,
        Pp,
        composition,
        cn_perfect,
        alpha,
    )

    coh_fcc["E"] = coh_fcc["E"].fillna(0)

    return {
        "rrange": rrange,
        "rhorange": rhorange,
        "rho": rho,
        "Fr": Fr,
        "Pp": Pp,
        "coh_fcc": coh_fcc,
        "element_energy_fcc": element_energy_fcc,
        "fcc_mean": float(coh_fcc["E"].iloc[0]),
        "fcc_std": float(coh_fcc["E"].iloc[1]),
    }


def calculate_vfe_statistics(
    unique_environments: Sequence[np.ndarray],
    frequencies: np.ndarray,
    lattice_parameter: float,
    composition: np.ndarray,
    alpha: np.ndarray,
    potential_data: Dict[str, Any],
) -> Dict[str, Any]:
    """Calculate analytical vacancy-formation-energy statistics."""
    number_of_sites = len(unique_environments)
    statistics = np.zeros((number_of_sites, 6), dtype=float)

    coh_fcc = potential_data["coh_fcc"]
    element_energy_fcc = potential_data["element_energy_fcc"]

    mean_vfe = 0.0
    variance_sum_1 = 0.0
    variance_sum_2 = 0.0

    for index, environment in enumerate(unique_environments):
        cn_local = copy.deepcopy(environment)
        cn_local[:, 0] *= lattice_parameter

        coh_site, _, element_energy_vfe = pot_defect.potential_stats(
            potential_data["rrange"],
            potential_data["rhorange"],
            potential_data["rho"],
            potential_data["Fr"],
            potential_data["Pp"],
            composition,
            cn_local,
            alpha,
        )

        coh_site["E"] = coh_site["E"].fillna(0)

        fcc_mean = float(coh_fcc["E"].iloc[0])
        fcc_std = float(coh_fcc["E"].iloc[1])
        site_mean = float(coh_site["E"].iloc[0])
        site_std = float(coh_site["E"].iloc[1])

        covariance = (
            np.sum(
                composition
                * element_energy_vfe
                * element_energy_fcc
            )
            - site_mean * fcc_mean
        )

        statistics[index, :] = [
            index,
            fcc_mean,
            fcc_std,
            site_mean,
            site_std,
            covariance,
        ]

        frequency = frequencies[index]

        mean_vfe += (
            site_mean - fcc_mean
        ) * frequency

        # Preserves the variance expression from the original code.
        sigma_excess = (
            site_std**2
            + fcc_std**2
            - 2.0 * site_std * fcc_std
        )

        variance_sum_1 += frequency * sigma_excess

        variance_sum_2 += (
            2.0
            * sigma_excess
            * frequency
            * (frequency - 1)
            / 2.0
        )

    std_vfe = np.sqrt(
        variance_sum_1 + variance_sum_2
    )

    statistics_df = pd.DataFrame(
        statistics,
        columns=[
            "VFE atom site",
            "Efcc average",
            "Efcc stdev",
            "Evfe average",
            "Evfe stdev",
            "covar",
        ],
    )

    std_vfe_check = np.sqrt(
        np.sum(
            frequencies.astype(float) ** 2
            * (
                statistics_df["Efcc stdev"].to_numpy() ** 2
                + statistics_df["Evfe stdev"].to_numpy() ** 2
                - 2.0
                * statistics_df["Evfe stdev"].to_numpy()
                * statistics_df["Efcc stdev"].to_numpy()
            )
        )
    )

    return {
        "mean_vfe": float(mean_vfe),
        "std_vfe": float(std_vfe),
        "std_vfe_check": float(std_vfe_check),
        "site_statistics": statistics_df,
    }


def calculate_ts_vme_statistics(
    unique_environments: Sequence[np.ndarray],
    frequencies: np.ndarray,
    lattice_parameter: float,
    composition: np.ndarray,
    alpha: np.ndarray,
    potential_data: Dict[str, Any],
    mean_vfe: float,
    std_vfe: float,
) -> Dict[str, Any]:
    """
    Calculate transition-state excess-energy and VME statistics.
    """
    number_of_sites = len(unique_environments)
    statistics = np.zeros((number_of_sites, 6), dtype=float)

    coh_fcc = potential_data["coh_fcc"]

    mean_ts = 0.0
    variance_sum_1 = 0.0
    variance_sum_2 = 0.0

    for index, environment in enumerate(unique_environments):
        cn_local = copy.deepcopy(environment)
        cn_local[:, 0] *= lattice_parameter

        coh_site, _, _ = pot_defect.potential_stats(
            potential_data["rrange"],
            potential_data["rhorange"],
            potential_data["rho"],
            potential_data["Fr"],
            potential_data["Pp"],
            composition,
            cn_local,
            alpha,
        )

        coh_site["E"] = coh_site["E"].fillna(0)

        fcc_mean = float(coh_fcc["E"].iloc[0])
        fcc_std = float(coh_fcc["E"].iloc[1])
        site_mean = float(coh_site["E"].iloc[0])
        site_std = float(coh_site["E"].iloc[1])

        statistics[index, :] = [
            index,
            fcc_mean,
            fcc_std,
            site_mean,
            site_std,
            0.0,
        ]

        frequency = frequencies[index]

        mean_ts += (
            site_mean - fcc_mean
        ) * frequency

        # Preserves the variance expression from the original code.
        sigma_excess = (
            site_std**2
            + fcc_std**2
            - 2.0 * site_std * fcc_std
        )

        variance_sum_1 += frequency * sigma_excess

        variance_sum_2 += (
            2.0
            * sigma_excess
            * frequency
            * (frequency - 1)
            / 2.0
        )

    std_ts = np.sqrt(
        variance_sum_1 + variance_sum_2
    )

    mean_vme = mean_ts - mean_vfe
    std_vme = np.sqrt(
        std_ts**2 + std_vfe**2
    )

    statistics_df = pd.DataFrame(
        statistics,
        columns=[
            "TS atom site",
            "Efcc average",
            "Efcc stdev",
            "ETS average",
            "ETS stdev",
            "covar",
        ],
    )

    return {
        "mean_ts": float(mean_ts),
        "std_ts": float(std_ts),
        "mean_vme": float(mean_vme),
        "std_vme": float(std_vme),
        "site_statistics": statistics_df,
    }


def run_vfe_vme_calculation(
    *,
    potential_file: PathLike,
    lattice_parameter: float,
    composition: Sequence[float],
    fcc_coordination_file: PathLike,
    vacancy_environment_file: PathLike,
    transition_state_environment_file: PathLike,
    use_alpha: bool = False,
    alpha_file: PathLike | None = None,
    environment_normalization_length: float = 3.49869654884664,
) -> Dict[str, Any]:
    """
    Run the complete FCC, VFE, and VME analytical calculation.

    The perfect-FCC coordination array is loaded from a pickle file.
    Therefore, ``rdf_coord.py`` and the perfect LAMMPS dump are not
    required during the calculation.
    """
    composition_array = np.asarray(
        composition,
        dtype=float,
    )

    if composition_array.ndim != 1:
        raise ValueError(
            "Composition must be one-dimensional."
        )

    if np.any(composition_array < 0):
        raise ValueError(
            "Composition fractions cannot be negative."
        )

    if not np.isclose(
        composition_array.sum(),
        1.0,
        atol=1.0e-6,
    ):
        raise ValueError(
            "Composition fractions must sum to 1. "
            f"Current sum = {composition_array.sum():.10f}"
        )

    cn_fcc = load_fcc_coordination(
        fcc_coordination_file
    )

    alpha = load_alpha_parameters(
        alpha_file=alpha_file,
        use_alpha=use_alpha,
        number_of_shells=cn_fcc.shape[0],
        number_of_elements=len(composition_array),
    )

    potential_data = calculate_fcc_reference(
        potential_file=potential_file,
        composition=composition_array,
        lattice_parameter=lattice_parameter,
        cn_fcc=cn_fcc,
        alpha=alpha,
    )

    vacancy_environments = (
        load_coordination_environments(
            pickle_file=vacancy_environment_file,
            normalization_length=(
                environment_normalization_length
            ),
        )
    )

    (
        unique_vacancy_environments,
        vacancy_frequencies,
    ) = get_unique_environments(
        vacancy_environments
    )

    vfe_results = calculate_vfe_statistics(
        unique_environments=(
            unique_vacancy_environments
        ),
        frequencies=vacancy_frequencies,
        lattice_parameter=lattice_parameter,
        composition=composition_array,
        alpha=alpha,
        potential_data=potential_data,
    )

    transition_state_environments = (
        load_coordination_environments(
            pickle_file=(
                transition_state_environment_file
            ),
            normalization_length=(
                environment_normalization_length
            ),
        )
    )

    (
        unique_ts_environments,
        ts_frequencies,
    ) = get_unique_environments(
        transition_state_environments
    )

    vme_results = calculate_ts_vme_statistics(
        unique_environments=unique_ts_environments,
        frequencies=ts_frequencies,
        lattice_parameter=lattice_parameter,
        composition=composition_array,
        alpha=alpha,
        potential_data=potential_data,
        mean_vfe=vfe_results["mean_vfe"],
        std_vfe=vfe_results["std_vfe"],
    )

    return {
        "E_FCC": potential_data["fcc_mean"],
        "Sig_FCC": potential_data["fcc_std"],
        "E_VFE": vfe_results["mean_vfe"],
        "Sig_VFE": vfe_results["std_vfe"],
        "Sig_VFE_check": (
            vfe_results["std_vfe_check"]
        ),
        "E_TS": vme_results["mean_ts"],
        "Sig_TS": vme_results["std_ts"],
        "E_VME": vme_results["mean_vme"],
        "Sig_VME": vme_results["std_vme"],
        "FCC_statistics": potential_data["coh_fcc"],
        "VFE_site_statistics": (
            vfe_results["site_statistics"]
        ),
        "TS_site_statistics": (
            vme_results["site_statistics"]
        ),
        "alpha": alpha,
        "cn_FCC": cn_fcc,
        "cn_VFE_unique": (
            unique_vacancy_environments
        ),
        "cn_VFE_freq": vacancy_frequencies,
        "cn_TS_unique": unique_ts_environments,
        "cn_TS_freq": ts_frequencies,
    }


def print_summary(
    results: Dict[str, Any],
    precision: int = 8,
) -> None:
    """
    Print FCC cohesive energy and defect-energy results.
    """
    print("=" * 52)
    print("ANALYTICAL ENERGY SUMMARY")
    print("=" * 52)

    print(
        f"FCC cohesive energy average = "
        f"{results['E_FCC']:.{precision}f} eV"
    )
    print(
        f"FCC cohesive energy stdev   = "
        f"{results['Sig_FCC']:.{precision}f} eV"
    )

    print("-" * 52)

    print(
        f"VFE average                 = "
        f"{results['E_VFE']:.{precision}f} eV"
    )
    print(
        f"VFE stdev                   = "
        f"{results['Sig_VFE']:.{precision}f} eV"
    )

    print("-" * 52)

    print(
        f"TS excess energy average    = "
        f"{results['E_TS']:.{precision}f} eV"
    )
    print(
        f"TS excess energy stdev      = "
        f"{results['Sig_TS']:.{precision}f} eV"
    )

    print("-" * 52)

    print(
        f"VME average                 = "
        f"{results['E_VME']:.{precision}f} eV"
    )
    print(
        f"VME stdev                   = "
        f"{results['Sig_VME']:.{precision}f} eV"
    )

    print("=" * 52)

# =========================================================
# RANDOM-ALLOY GPFE FUNCTIONS
# =========================================================

GPFE_FAULT_TYPES = (
    "USF",
    "ISF",
    "UTF1",
    "ESF",
    "UTF2",
    "TF",
)


def calculate_random_gpfe(
    *,
    potential_file: PathLike,
    lattice_parameter: float,
    cutoff_radius: float,
    composition: Sequence[float],
    fcc_coordination_file: PathLike,
    fault_types: Sequence[str] = GPFE_FAULT_TYPES,
) -> Dict[str, Any]:
    """
    Calculate the random-alloy generalized planar fault energies.

    Notes
    -----
    This function intentionally does not accept Warren-Cowley parameters.
    It uses the original random-alloy ``Potential_GPFE.py`` implementation.

    The normalized FCC coordination array is loaded from ``cn_FCC.pkl``.
    Faulted coordination arrays are generated using
    ``rdf_coord.rdf_coord_fault``.
    """
    composition_array = np.asarray(composition, dtype=float)

    if composition_array.ndim != 1:
        raise ValueError("Composition must be one-dimensional.")

    if np.any(composition_array < 0):
        raise ValueError("Composition fractions cannot be negative.")

    if not np.isclose(composition_array.sum(), 1.0, atol=1.0e-6):
        raise ValueError(
            "Composition fractions must sum to 1. "
            f"Current sum = {composition_array.sum():.10f}"
        )

    cn_fcc_normalized = load_fcc_coordination(
        fcc_coordination_file
    )

    # Convert normalized shell distances to angstrom for the GPFE routines.
    cn_fcc = copy.deepcopy(cn_fcc_normalized)
    cn_fcc[:, 0] *= lattice_parameter

    rrange, rhorange, rho, Fr, Pp = (
        pot_gpfe.potential_read(str(potential_file))
    )

    form_E_fcc, covar_fcc, E_element_fcc = (
        pot_gpfe.potential_stats(
            rrange,
            rhorange,
            rho,
            Fr,
            Pp,
            composition_array,
            cn_fcc,
        )
    )

    fault_results: Dict[str, pd.DataFrame] = {}
    fault_covariances: Dict[str, float] = {}
    fault_coordination: Dict[str, Any] = {}

    valid_fault_types = set(GPFE_FAULT_TYPES)

    for fault_type in fault_types:
        if fault_type not in valid_fault_types:
            raise ValueError(
                f"Unknown GPFE fault type '{fault_type}'. "
                f"Choose from {GPFE_FAULT_TYPES}."
            )

        _, cn_fault = rc.rdf_coord_fault(
            lattice_parameter,
            cutoff_radius,
            cn_fcc,
            fault_type,
        )

        form_E_fault, covar_fault = (
            pot_gpfe.potential_stats_fault(
                rrange,
                rhorange,
                rho,
                Fr,
                Pp,
                composition_array,
                cn_fcc,
                cn_fault,
                form_E_fcc,
                E_element_fcc,
                fault_type,
            )
        )

        fault_results[fault_type] = form_E_fault
        fault_covariances[fault_type] = float(covar_fault)
        fault_coordination[fault_type] = cn_fault

    summary_rows = []

    for fault_type in fault_types:
        table = fault_results[fault_type]
        column_name = f"E_{fault_type}"

        summary_rows.append(
            {
                "Fault": fault_type,
                "Mean (eV/atom)": float(
                    table.loc["Mean", column_name]
                ),
                "Std (eV/atom)": float(
                    table.loc["Std", column_name]
                ),
            }
        )

    summary_table = pd.DataFrame(summary_rows)

    return {
        "FCC_statistics_GPFE": form_E_fcc,
        "FCC_covariance_GPFE": float(covar_fcc),
        "FCC_element_energies_GPFE": E_element_fcc,
        "GPFE_results": fault_results,
        "GPFE_covariances": fault_covariances,
        "GPFE_summary": summary_table,
        "GPFE_fault_coordination": fault_coordination,
        "cn_FCC_GPFE": cn_fcc,
    }


def print_gpfe_summary(
    results: Dict[str, Any],
    precision: int = 8,
) -> None:
    """Print random-alloy GPFE means and standard deviations."""
    print("=" * 58)
    print("RANDOM-ALLOY GPFE SUMMARY")
    print("=" * 58)

    fcc_table = results["FCC_statistics_GPFE"]

    print(
        f"FCC cohesive energy average = "
        f"{float(fcc_table.loc['Mean', 'E']):.{precision}f} eV/atom"
    )
    print(
        f"FCC cohesive energy stdev   = "
        f"{float(fcc_table.loc['Std', 'E']):.{precision}f} eV/atom"
    )

    print("-" * 58)

    for _, row in results["GPFE_summary"].iterrows():
        print(
            f"{row['Fault']:<5s} mean = "
            f"{row['Mean (eV/atom)']:.{precision}f} eV/atom, "
            f"stdev = {row['Std (eV/atom)']:.{precision}f} eV/atom"
        )

    print("=" * 58)
