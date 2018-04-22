# MECPy

Python utilities for performing MECP (Minimum Energy Crossing Point) with Gaussian.

MECPy provides a Python wrapper around the [original J. N. Harvey's](https://link.springer.com/article/10.1007/s002140050309) Fortran code that greatly simplifies the setup and includes small improvements gathered from other MECP variants ([sobMECP](http://sobereva.com/286)).

_For the sake of completeness, I have also included MECPro: a MIT-licensed, Python-only fork of MECP, originally developed by [J. Snyder and D. H. Ess at BYU](http://jur.byu.edu/?p=22227), updated so it is pip-installable and unit-testable._

## Features

- You don't have to modify the source code anymore. All the setup steps are performed automatically. For example, the number of atoms is inferred from the geometry input.
- Geometry can be specified with atomic symbols or numbers. The script will take care of the rest.
- Self-contained Python executable with no 3rd party dependencies: you can move it around!
- It still uses the original MECP code behind the scenes, so you can trust the results.
- Optimization trajectory is written for every step.
- Several energy parsers included: `dft`, `pm2`, `ci`, `td`.
- Unit-tests (needs `pytest`).

## Requirements

- Python 2.7/3.4+
- `gfortran` or equivalent compiler

## Why?!

The original MECP code requires a rather intricate setup process. For example, the user is expected to manually modify the code and recompile the Fortran binaries for each calculation. My opinion is that it should be easier and that's why I automated the manual steps in a script called `mecpsh`.

During the development of MECPsh (the Python wrapper around MECP), I found about MECPro (a Python-only fork), which I have included here for completeness, in a pip-installable form. The authors claimed reproducibility of the calculations, but I could not locate further evidence, so I also included some unit-tests extracted from ([sobMECP](http://sobereva.com/286)) [WIP].

## Usage

The `mecpy` package includes a self-contained, pure-Python module called `mecpsh` that can be called directly with `python mecpsh.py`. If the package is installed with `pip`, the executable `mecpsh` will also be available. Either way, the `-h` flag will show the usage:

```
# Direct execution
python mecpsh.py -h
# Only available if installed with pip
mecpsh -h
```

All the configuration keys have default values, so, if all the requested files are named like that, all you have to do is run `python mecpsh.py` in the corresponding directory. If custom values are preferred, you can specify them with the appropriate flags OR through a configuration file.

```
python mecpsh.py --geom initial_geometry --FC ifort
```


### Required files

- `Input_Header_A`: File containing the top part of system configuration with multiplicity A. Name can be changed with `--a_header` key.
- `Input_Header_B`: File containing the top part of system configuration with multiplicity B. Name can be changed with `--b_header` key.
- `geom`: File containing the starting system geometry (element symbols will be converted in atomic numbers automatically). Name can be changed with `--geom`` key.
- `footer` (optional): File containing the bottom part of system configuration (basis sets, etc). Name can be changed with `--footer` flag.

_More details can be found in the command line help message with `python mecpsh.py -h` or `mecpsh -h`._


### Output

All the output files will be written to the working directory and stored in the `JOBS` folder, except the one corresponding to the last step. That folder will also contain a `trajectory.xyz` step with all the intermediate geometries. As in the original MECP, `ReportFile` will list more information about the optimization process
