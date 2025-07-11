# Johnson–Cook VUMAT with Damage and Taylor–Quinney Heating
[![DOI](https://zenodo.org/badge/719144148.svg)](https://doi.org/10.5281/zenodo.15839293) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This repository contains a user-defined material subroutine (VUMAT) implementing the Johnson–Cook constitutive model for finite deformation plasticity with strain rate and temperature dependence, as well as isotropic damage evolution. The code is developed for use in `Abaqus/Explicit` and written in double precision Fortran.

The implementation is fully documented and based on a return mapping algorithm using the bisection method for stable plastic correction. Damage initiation is based on the Johnson–Cook criterion, while damage evolution follows a displacement-based linear softening law (Displacement at Failure).

## Documentation

A detailed derivation of the governing equations and numerical implementation is provided in the following PDF:

**[johnson-cook-vumat-doc.pdf](./johnson-cook-vumat-doc.pdf)**

## Validation

The `validation/` folder contains input files and plots for comparing this VUMAT implementation against Abaqus' built-in Johnson–Cook model. The test consists of a single-element cube made of Ti-6Al-4V undergoing uniaxial deformation.

### Input Files

- `jc_ti64_abq_builtin.inp` — Abaqus built-in Johnson–Cook model
- `jc_ti64_vumat.inp` — Custom Johnson–Cook VUMAT

### Comparison Plots

The plots below compare the stress–strain response and damage evolution between the two models. The agreement confirms the accuracy of the implementation.

<p align="center">
  <img src="validation/stress_strain_comparison.png" width="500"/>
</p>

<p align="center">
  <img src="validation/damage_criterion_comparison.png" width="500"/>
</p>

### Citations

Material parameters used in the validation test are referenced in `validation/citations.md`.

## Cite this repository

If you use this code in your work, please cite it using the following Zenodo record:

**DOI:** [10.5281/zenodo.15839293](https://doi.org/10.5281/zenodo.15839293)

BibTeX:
```
@software{arcidiacono2025johnsoncook,
  author       = {Arcidiacono, Mauro F. and Rahimi, Salaheddin},
  title        = {Johnson–Cook VUMAT with Damage and Thermal Effects},
  version      = {v1.0.0},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.15839293},
  url          = {https://doi.org/10.5281/zenodo.15839293}
}
```
