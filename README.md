# Johnson–Cook VUMAT with Damage and Thermal Effects

This repository contains a user-defined material subroutine (VUMAT) implementing the Johnson–Cook constitutive model for finite deformation plasticity with strain rate and temperature dependence, as well as isotropic damage evolution. The code is developed for use in `Abaqus/Explicit` and written in double precision Fortran.

The implementation is fully documented and based on a return mapping algorithm using the bisection method for stable plastic correction.

## Documentation

A detailed derivation of the governing equations and numerical implementation is provided in the following PDF:

**[johnson-cook-vumat-derivation.pdf](./johnson-cook-vumat-derivation.pdf)**
