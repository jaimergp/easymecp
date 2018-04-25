# MECPy

Python utilities for performing MECP (Minimum Energy Crossing Point) with Gaussian.

MECPy provides a Python wrapper around the [original J. N. Harvey's](https://link.springer.com/article/10.1007/s002140050309) Fortran code that greatly simplifies the setup and includes small improvements gathered from other MECP variants ([sobMECP](http://sobereva.com/286)).

_For the sake of completeness, I have also included MECPro: a MIT-licensed, Python-only fork of MECP, originally developed by [J. Snyder and D. H. Ess at BYU](http://jur.byu.edu/?p=22227), updated so it is pip-installable and unit-testable._

## Features

- You don't have to modify the source code anymore. All the setup steps are performed automatically. For example, the number of atoms is inferred from the geometry input.
- Geometry can be specified with atomic symbols or numbers. The script will take care of the rest.
- Convergence criteria can be easily modified.
- No hardcoded values: use another Gaussian version, Fortran compiler, flags...
- Self-contained Python executable with no 3rd party dependencies: you can move it around!
- It still uses the original MECP code behind the scenes, so you can trust the results.
- Optimization trajectory is written for every step.
- Several energy parsers included: `dft`, `pm2`, `ci`, `td`. You can include a new one in a separate Python file if needed.
- Unit-tests (needs `pytest`).

## Why?!

The original MECP code requires a rather intricate setup process. For example, the user is expected to manually modify the code and recompile the Fortran binaries for each calculation. My opinion is that it should be easier and that's why I automated the manual steps in a script called `mecpy`.

During the development of MECPy (the Python wrapper around MECP), I found about MECPro (a Python-only fork), which I have included here for completeness, in a pip-installable form. The authors claimed reproducibility of the calculations, but I could not locate further evidence, so I also included some unit-tests extracted from ([sobMECP](http://sobereva.com/286)) [WIP].

## Installation

If you just want to run a few calculations, simply grab the latest `mecpy.py` script from [here](https://github.com/jaimergp/mecpy/blob/master/mecpy/mecpy.py). 

If you'd like a permanent installation, use `pip` like this:

```
pip install https://github.com/jaimergp/mecpy/archive/master.zip
```

### Requirements

- Python 2.7/3.4+
- `gfortran` or equivalent compiler
- `pytest` (only for running tests)


## Usage

The `mecpy` package includes a self-contained, pure-Python module called `mecpy` that can be called directly with `python mecpy.py`. If the package is installed with `pip`, the executable `mecpy` will also be available. Either way, the `-h` flag will show the usage:

```
# Direct execution
python mecpy.py -h
# Only available if installed with pip
mecpy -h
```

All the configuration keys have default values, so, if all the requested files are named like that, all you have to do is run `python mecpy.py` in the corresponding directory. If custom values are preferred, you can specify them with the appropriate flags OR through a configuration file.

```
python mecpy.py --geom initial_geometry --FC ifort
```

### Required files

- `Input_Header_A`: File containing the top part of system configuration with multiplicity A. Name can be changed with `--a_header` key.
- `Input_Header_B`: File containing the top part of system configuration with multiplicity B. Name can be changed with `--b_header` key.
- `geom`: File containing the starting system geometry (element symbols will be converted in atomic numbers automatically). Name can be changed with `--geom`` key.
- `footer` (optional): File containing the bottom part of system configuration (basis sets, etc). Name can be changed with `--footer` flag.

_More details can be found in the command line help message with `python mecpy.py -h` or `mecpy -h`._

__TIP__: Use `--gaussian_exe` key to specify the Gaussian executable (version) to use: `g09` or `g16`. Others might work as well. If not specified, `mecpy` will use the default value: `g16` if Gaussian 16 is present in `$PATH`; `g09` otherwise. You can also specify absolute paths here, if needed.


### Output

All the output files will be written to the working directory and stored in the `JOBS` folder, except the one corresponding to the last step. That folder will also contain a `trajectory.xyz` step with all the intermediate geometries. As in the original MECP, `ReportFile` will list more information about the optimization process


### Examples

Run a MECP calculation providing all files with default names in the same directory as `mecpy.py`:

```
python mecpy.py
```

Increase the maximum allowed steps:

```
python mecpy.py --max_steps 100
```

Less tight TDE convergence criteria (you must use Fortran-style double-precision floats)

```
python mecpy.py --TDE 5.d-4
```

Provide a different energy parser (example file in `tests/data/energy.py`; must include a `parse_energy` toplevel function)

```
python mecpy.py --energy_parser energy.py
```

### Restarting jobs

If the calculation does not converge before reaching `max_steps`, you might want to extend it. It's simple: take the latest geometry you like (last one is always in the file `geom`, but you can create another if you want), and relaunch `mecpy --geom <your_geometry>`. Since the `JOBS` folder is automatically renamed to `JOBS1`, `JOBS2`, etc, you won't lose your files.
