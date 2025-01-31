# WTB

This is the improved version of [WanTiBexos](https://github.com/ac-dias/wantibexos) code,

see [DOI](https://doi.org/10.1016/j.cpc.2022.108636)

## Version 1.0
What's new:

- the MKL solver for Bethe-Salpeter equation was replaced by MPI+GPU paralel solver from the Scalable Library for Eigenvalue Problem Computations (SLEPc)

The calculation time for Example-3 is reduced from 680 s using MKL solver to 3 s using SLEPc solver on 64 cores.

How to install the code:

  
