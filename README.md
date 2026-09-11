# Solid-Solution Energy Statistics

Analytical calculation of the statistical distributions of energies in concentrated solid solutions (CSSs) described by embedded-atom-method (EAM) potentials. Given an EAM potential, a composition, a lattice parameter, and (optionally) Warren-Cowley short-range-order parameters, the code returns the mean and standard deviation of:

* the per-atom cohesive energy of the perfect FCC crystal, for random and short-range ordered (SRO) alloys;
* the vacancy formation energy (VFE), for random and SRO alloys;
* the transition-state (TS) excess energy and the vacancy migration energy (VME), for random and SRO alloys;
* the generalized planar fault energies (GPFE: USF, ISF, UTF1, ESF, UTF2, TF), for random alloys.

No simulations are run. Each calculation takes about a second on a laptop and uses only the coefficients of the interatomic potential, the alloy composition, and the coordination structure factors of the perfect and defective lattice.

<p align="center">
  <img src="Images/Flowchart.png" alt="Program workflow" width="750">
</p>

---

## Citation

The analytical framework is described in two papers. Please cite both when using the code:

> R. Jagatramka, C. Wang, and M. Daly, "An analytical method to quantify the statistics of energy landscapes in random solid solutions," *Computational Materials Science* **214** (2022) 111763. [doi:10.1016/j.commatsci.2022.111763](https://doi.org/10.1016/j.commatsci.2022.111763)
> (cohesive-energy statistics and generalized planar fault energies for random alloys)

> A. Baski, R. Jagatramka, and M. Daly, "A mechanistic model for vacancy energetics in concentrated solid solutions with short-range order," in preparation (2026).
> (short-range order, vacancy formation energy, and vacancy migration energy)

A machine-readable [`CITATION.cff`](CITATION.cff) is included; the second entry will be updated on publication. The code is distributed under the [GPL-3.0 license](LICENSE).

---

## Requirements and installation

* Python 3.8 or later
* NumPy and pandas (`pip install -r requirements.txt`)

```bash
git clone https://github.com/dalymatt/Solid-Solution-PE-Statistics
cd Solid-Solution-PE-Statistics
pip install -r requirements.txt
python example_defect_energy.py
```

The last command runs the shipped example (Fe73Ni8Cr19, SRO mode) and should finish in a few seconds. Expected output is given under [Reference values](#reference-values).

---

## Running a calculation

All inputs live in the "MATERIAL-SPECIFIC INPUTS" block of `example_defect_energy.py`. Edit that block and run the script; nothing else needs to change.

```python
POTENTIAL_FILE    = ROOT / "potentials" / "FeNiCr.eam.alloy"
LATTICE_PARAMETER = 3.51036                        # angstrom
CUTOFF_RADIUS     = 5.6                            # angstrom, from the potential file
COMPOSITION       = np.array([0.73, 0.08, 0.19])   # mole fractions, order as in the potential
MODE              = "Random"                       # "Random" or "SRO"
ALPHA_FILE        = ROOT / "data" / "sro" / "304SS" / "alpha_Fe73Ni8Cr19_alpha_p0p05.npy"   # SRO only
```

| Input | Notes |
|---|---|
| `POTENTIAL_FILE` | An EAM/alloy potential in LAMMPS `setfl` format, stored in `potentials/`. Three are shipped: `FeNiCr.eam.alloy` (Bonny et al.), `FeNiCrCoCu-with-ZBL.eam.alloy` (Deluigi et al.), `NiCo-lammps-2014.alloy` (Béland et al.). Others can be obtained from the [NIST Interatomic Potentials Repository](https://www.ctcms.nist.gov/potentials/). |
| `LATTICE_PARAMETER` | FCC lattice parameter of the alloy at the composition and SRO state of interest, in Å. The values used in the papers are listed in comments in the example script. |
| `CUTOFF_RADIUS` | Cutoff of the potential, in Å (5.6 for FeNiCr, 5.80375 for FeNiCrCoCu, 6.5 for NiCo). It sets how many FCC neighbour shells are included (5 for all shipped cases). |
| `COMPOSITION` | Mole fractions summing to 1, in the element order of the potential file. A five-element potential can describe a ternary by setting the unused fractions to 0. |
| `MODE` | `"Random"` computes cohesive, VFE/VME, and GPFE statistics. `"SRO"` computes cohesive and VFE/VME statistics with Warren-Cowley parameters; GPFE is skipped. |
| `ALPHA_FILE` | Warren-Cowley parameter array, used only when `MODE = "SRO"`. Format and shipped files are described under [Warren-Cowley parameter files](#warren-cowley-parameter-files). |

The script prints three blocks: the FCC cohesive-energy table (mean and standard deviation of the electron density `rho`, embedding energy `F`, pair-interaction energy `Pp`, and cohesive energy `E`, in eV/atom), the defect-energy summary (VFE, TS excess energy, VME, in eV), and, in Random mode, the GPFE tables in eV/atom and mJ/m². All results are also returned in the `results` dictionary for use in your own scripts (see the end of the example file).

---

## Reference values

Use these to confirm an installation. Lattice parameters are those used in the papers; cutoffs are 5.6 Å for the FeNiCr potential and 5.80375 Å for the FeNiCrCoCu potential.

| System | a (Å) | Mode | VFE (eV) | VME (eV) |
|---|---|---|---|---|
| Fe73Ni8Cr19 | 3.51036 | Random | 2.021 ± 0.032 | 0.951 ± 0.055 |
| Ni33Cr33Co34 | 3.53073 | Random | 1.455 ± 0.199 | 1.106 ± 0.211 |
| Ni25Cr25Co25Cu25 | 3.546 | Random | 1.321 ± 0.176 | 0.980 ± 0.219 |
| Fe20Ni20Cr20Co20Cu20 | 3.54939 | Random | 1.426 ± 0.225 | 1.041 ± 0.257 |
| Fe73Ni8Cr19 | 3.51036 | SRO, `alpha_p0p05` | 2.019 ± 0.032 | 0.944 ± 0.057 |

Full output of the shipped example (`MODE = "SRO"`):

```text
FCC COHESIVE-ENERGY STATISTICS: SRO SYSTEM
           rho         F        Pp         E
Mean  1.012853 -2.106122 -2.047309 -4.153431
Std   0.013594  0.329822  0.068708  0.282410

DEFECT ENERGY SUMMARY
VFE average                 = 2.01866673 eV
VFE stdev                   = 0.03190829 eV
TS excess energy average    = 2.96254142 eV
TS excess energy stdev      = 0.04746717 eV
VME average                 = 0.94387469 eV
VME stdev                   = 0.05719503 eV
```

The same alloy with `MODE = "Random"`:

```text
FCC COHESIVE-ENERGY STATISTICS: RANDOM SYSTEM
           rho         F        Pp         E
Mean  1.012853 -2.106131 -2.048497 -4.154628
Std   0.013401  0.329808  0.066111  0.284211

DEFECT ENERGY SUMMARY
VFE average                 = 2.02080912 eV
VFE stdev                   = 0.03242783 eV
TS excess energy average    = 2.97176968 eV
TS excess energy stdev      = 0.04437746 eV
VME average                 = 0.95096055 eV
VME stdev                   = 0.05496292 eV

GENERALIZED PLANAR FAULT ENERGY STATISTICS
         E_USF     E_ISF    E_UTF1     E_ESF    E_UTF2      E_TF
Mean  0.125037  0.005899  0.127986  0.005899  0.127986  0.005899
Std   0.062559  0.054747  0.062409  0.054942  0.062193  0.054747

GENERALIZED PLANAR FAULT ENERGY STATISTICS (mJ/m^2)
           E_USF       E_ISF      E_UTF1       E_ESF      E_UTF2        E_TF
Mean  375.442787   17.713483  384.299567   17.713478  384.299517   17.713483
Std   187.842268  164.387701  187.392144  164.973490  186.744363  164.387701
```

Two further cohesive-energy checks from the 2022 paper: Ni40Co60 (`NiCo-lammps-2014.alloy`, a = 3.512 Å, cutoff 6.5 Å) and Fe33Ni33Cr34 (`FeNiCr.eam.alloy`, a = 3.5225 Å, cutoff 5.6 Å).

<p align="center">
  <img src="Images/stats.png" alt="Ni40Co60 cohesive-energy statistics" width="500">
  <img src="Images/stats2.png" alt="Fe33Ni33Cr34 cohesive-energy statistics" width="500">
</p>

---

## Warren-Cowley parameter files

An SRO calculation reads a NumPy array of shape `(n_shells, n_elements, n_elements)` from `ALPHA_FILE`. Entry `[k, X, Y]` is the Warren-Cowley parameter α for finding species Y in neighbour shell k+1 around a central atom of species X, defined through the conditional probability p(Y|X) = c_Y (1 − α_XY). Element order follows the potential file; shell order follows the FCC shells inside `CUTOFF_RADIUS`, so `n_shells` must equal the number of shells the cutoff admits (5 for all shipped cases; a larger cutoff needs a larger array).

On loading, the code enforces the sum rule Σ_Y c_Y α_XY = 0 in every shell (tolerance 10⁻⁶) and requires p(Y|X) ≥ 0. A file that violates either condition raises an error naming the shell and pair. Because of the sum rule, the homoatomic entries are fixed by the heteroatomic ones; for identical heteroatomic values α this gives α_XX = −α (1 − c_X) / c_X, and a negative α is only possible if c_X ≥ −α / (1 − α) for every species.

Shell-resolved values are applied to the vacancy and saddle-point environments by distance: each coordination peak takes the α of the nearest ideal FCC shell.

Shipped files, in `data/sro/304SS/` (Fe73Ni8Cr19) and `data/sro/Quinary/` (Fe20Ni20Cr20Co20Cu20), were measured from Monte Carlo swapped structures targeting a first-shell heteroatomic value of −0.05, −0.025, 0, +0.025, +0.05, +0.10, and +0.15 (filename suffix `m0p05` … `p0p15`; `m` negative, `p` positive). Outer-shell values are emergent rather than targeted and can be large for dilute species. The Fe73Ni8Cr19 `p0p00` file is a sampled random structure and carries small non-zero outer-shell values; the Quinary `p0p00` file is exactly zero. The Monte Carlo code that produced these files is not part of this repository.

---

## Defect environment files

`data/coordination/cn_vac.pkl` and `cn_TS.pkl` hold, for every atom within the cutoff of a monovacancy (78 atoms) and of a migrating atom at the saddle point (456 atoms), the list of neighbour distances and multiplicities. They were extracted from relaxed pure-Fe molecular-statics and climbing-image nudged-elastic-band structures under the Bonny et al. FeNiCr potential, are stored normalized by that structure's lattice parameter (3.49869654884664 Å), and are rescaled to `LATTICE_PARAMETER` at run time. The same environments are used for every alloy, as described in Section 3.2 of Baski et al.

The perfect-FCC reference is not read from a file. It is generated from the ideal shell ratios (√½, 1, √(3/2), √2, √(5/2), …) and `LATTICE_PARAMETER`, so the reference and the cohesive-energy calculation always use the same structure factors. Because VFE and VME are sums of per-site excess energies over 78 and 456 sites, any inconsistency between the reference and the environment distances is multiplied by those numbers; environment files regenerated from a differently prepared structure will shift the means at the level of a few hundredths of an eV. The scripts that generated the environment files are not part of this repository.

`data/coordination/cn_USF.pkl`, `cn_ISF.pkl`, `cn_UTF1.pkl`, `cn_ESF.pkl`, `cn_UTF2.pkl`, and `cn_TF.pkl` hold the per-layer structure factors of the six planar faults used by the GPFE calculation (see the 2022 paper, Fig. 4 and Table 1).

---

## Generalized planar fault energies

In Random mode the code evaluates the six critical energies of the GPFE landscape: unstable stacking fault (USF), intrinsic stacking fault (ISF), first unstable twinning fault (UTF1), extrinsic stacking fault (ESF), second unstable twinning fault (UTF2), and twin fault (TF). Each is the excess energy per atom of the non-FCC {111} layers around the fault, summed over those layers, so the first table is in eV/atom. The second table converts to a planar fault energy in mJ/m² using the {111} areal density ρ₁₁₁ = 4 / (√3 a²) and 1 eV/Å² = 16021.77 mJ/m². Means and standard deviations follow Eqs. 11 to 16 of the 2022 paper.

Two reference calculations from the 2022 paper, in eV/atom:

<p align="center">
  <img src="Images/nicofaultedstate.PNG" alt="Ni40Co60 GPFE statistics" width="450">
  <img src="Images/frnicrfaultedstate.PNG" alt="Fe33Ni33Cr34 GPFE statistics" width="450">
</p>

Left: Ni40Co60, a = 3.512 Å, cutoff 6.5 Å. Right: Fe33Ni33Cr34, a = 3.5225 Å, cutoff 5.6 Å.

GPFE statistics for SRO alloys are not implemented; `MODE = "SRO"` skips this block.

---

## Repository layout

```text
example_defect_energy.py     the input script you edit and run
src/energy_workflow.py       orchestrates the cohesive, VFE/VME, and GPFE calculations
src/Potential.py             EAM setfl reader; cohesive-energy statistics with SRO
src/defect_energy_functions.py   VFE and VME statistics from the defect environments
src/Potential_GPFE.py        random-alloy cohesive and GPFE statistics (2022 paper)
src/rdf_coord.py             perfect-FCC structure factors; loader for faulted-layer structure factors
src/README_EQUATION_MAP.md   which equations each module implements
potentials/                  three EAM setfl files
data/coordination/           defect and planar-fault environment files
data/sro/                    Warren-Cowley parameter files for the two sweep systems
papers/                      the 2022 paper
```

`rdf_coord.py` also provides BCC and HCP shell tables (`rdf_coord_bcc`, `rdf_coord_hcp`) for use outside the workflow; the workflow itself is FCC only.

---

## Authors

Akash Baski, Ritesh Jagatramka, Chu Wang (now at Nissan), Ariana Sofia Del Valle, Amir Shirsalimian, and Matthew Daly. Advanced Materials and Microstructures Laboratory, University of Illinois Chicago. Corresponding author: mattdaly@uic.edu.

This work was supported by the National Science Foundation under Grant No. DMR-2144451.
