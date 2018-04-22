# MECPy

Python utilities for performing MECP (Minimum Energy Crossing Point) with Gaussian.

For the sake of completeness, this project includes several ways to perform MECP with Gaussian.

- MECP: The [original J. N. Harvey's](https://link.springer.com/article/10.1007/s002140050309) Fortran code.
- MECPsh: A Python wrapper around Fortran MECP that greatly simplifies the setup and includes small improvements gathered from other MECP variants ([sobMECP](http://sobereva.com/286)).
- MECPro: A Python-only fork of MECP, originally developed by [J. Snyder and D. H. Ess at BYU](http://jur.byu.edu/?p=22227), updated so it is pip-installable and unit-testable.
- MECPy: My own Python-only fork of MECP [WIP].

## Why?!

The original MECP code requires a rather intricate setup process. For example, the user is expected to manually modify the code and recompile the Fortran binaries for each calculation. My opinion is that it should be easier and that's why I coded the manual steps in a script called `mecpsh`.

During the development of MECPsh (the Python wrapper around MECP), I found about MECPro (a Python-only fork), which I have included here for completeness, in a pip-installable form. The authors claimed reproducibility of the calculations, but I could not locate further evidence, so I also included some unit-tests extracted from ([sobMECP](http://sobereva.com/286).

For educational purposes, I might end up coding a Python-only fork of my own, which I will call MECPy.