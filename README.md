# Solid-Solution Potential-Energy Statistics

Analytical calculations of cohesive-energy and defect-energy statistics in concentrated solid solutions using coordination-reparameterized EAM relations.

The repository supports:

- **Random systems:** FCC cohesive energy, vacancy formation energy (VFE), vacancy migration energy (VME), and generalized planar fault energies (GPFEs).
- **Short-range ordered systems:** FCC cohesive energy, VFE, and VME using Warren–Cowley short-range-order parameters.

## Quick start

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Open the root-level file:

```text
example_defect_energy.py
```

Select the calculation mode:

```python
MODE = "Random"
```

or

```python
MODE = "SRO"
```

Then run:

```bash
python example_defect_energy.py
```

The script uses paths relative to the repository location, so it can be launched from any working directory.

## Main inputs

Edit the following values in `example_defect_energy.py`:

```python
POTENTIAL_FILE
LATTICE_PARAMETER
CUTOFF_RADIUS
COMPOSITION
MODE
ALPHA_FILE
ENVIRONMENT_NORMALIZATION_LENGTH
```

The composition must follow the element order stored in the selected EAM/alloy potential file and must sum to one.

## Repository layout

```text
Solid-Solution-PE-Statistics-Organized/
├── example_defect_energy.py      # Main user-facing calculation file
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── src/                          # Analytical implementation
├── potentials/                   # EAM/alloy potential files
├── structures/                   # LAMMPS atomic structures
├── data/
│   ├── coordination/             # FCC, vacancy, TS, and fault coordination data
│   └── sro/                      # Warren–Cowley alpha arrays
├── examples/                     # Additional and legacy examples
├── papers/                       # Associated papers/manuscripts
├── images/                       # README and workflow figures
└── docs/                         # Additional documentation
```

## Calculation behavior

| Mode | FCC | VFE | VME | GPFE |
|---|---:|---:|---:|---:|
| Random | Yes | Yes | Yes | Yes |
| SRO | Yes | Yes | Yes | Not currently included |

When `MODE = "SRO"`, the main script automatically skips the GPFE calculation.

## Notes on supplied data

- `data/coordination/cn_FCC.pkl` contains the normalized perfect-FCC coordination structure.
- `data/coordination/cn_vac.pkl` contains vacancy-affected environments.
- `data/coordination/TS_rel.pkl` contains migration transition-state environments.
- `data/sro/*.npy` contains example Warren–Cowley SRO arrays.
- Potential files are stored separately under `potentials/`.
- LAMMPS configurations are stored separately under `structures/` and are not interatomic potentials.

## References

The random-solid-solution framework originates from:

R. Jagatramka, C. Wang, and M. Daly, “An analytical method to quantify the statistics of energy landscapes in random solid solutions,” *Computational Materials Science* 214 (2022) 111763.

The vacancy and SRO extension is documented in the manuscript included under `papers/`.

## License

See `LICENSE`.
