# Commented analytical-energy Python modules

The calculations are unchanged; comments and docstrings were added to identify
the physical quantities and manuscript equations.

## Equation map

- `Potential.py`: vacancy manuscript Eqs. (1)-(9); random-alloy counterparts in Jagatramka et al. (2022), Eqs. (1)-(10).
- `defect_energy_functions.py`: average VFE/VME, vacancy manuscript Eqs. (10)-(12); variances, Eqs. (13)-(14); GPFE, Jagatramka et al. (2022), Eqs. (11)-(16).
- `Potential_GPFE.py`: Jagatramka et al. (2022), cohesive Eqs. (1)-(10), GPFE Eqs. (11)-(16), and Appendices B-C.
- `energy_workflow.py`: orchestration layer linking the calculations above.

## References

1. A. Baski et al., *A mechanistic model for vacancy energetics in concentrated solid solutions with short-range order*.
2. R. Jagatramka, C. Wang, and M. Daly, *Computational Materials Science* 214 (2022) 111763. DOI: 10.1016/j.commatsci.2022.111763.
