# Input guide

The primary user interface is `example_defect_energy.py` in the repository root.

## Random calculation

```python
MODE = "Random"
```

No alpha file is used. FCC, VFE, VME, and GPFE statistics are calculated.

## SRO calculation

```python
MODE = "SRO"
ALPHA_FILE = ROOT / "data" / "sro" / "alpha_FeNiCr_SS_minus_point05.npy"
```

FCC, VFE, and VME statistics are calculated using the supplied Warren–Cowley array. GPFE is skipped.
