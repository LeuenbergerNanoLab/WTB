# WTB

This is the improved version of [WanTiBexos](https://github.com/ac-dias/wantibexos) code,

see A.C. Dias, J.F.R.V. Silveira, F. Qu "WanTiBEXOS: A Wannier based Tight Binding code for electronic band structure, excitonic and optoelectronic properties of solids" Comp. Phys. Comm., 285, 108636 (2022) [DOI](https://doi.org/10.1016/j.cpc.2022.108636)

## Version 1.0
What's new:

- the MKL solver for Bethe-Salpeter equation was replaced by MPI+GPU paralel solver from the [Scalable Library for Eigenvalue Problem Computations (SLEPc)](https://slepc.upv.es)

The calculation time for [Example-3](examples/example-3) is reduced from 680 s using MKL solver to 3 s using SLEPc solver on 64 cores.

How to install the code:

  
