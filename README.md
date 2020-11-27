# EasyMECP

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4293421.svg)](https://doi.org/10.5281/zenodo.4293421)


Self-contained Python script for performing MECP (Minimum Energy Crossing Point) with Gaussian.

easyMECP provides a single-file, self-contained Python wrapper around the [original J. N. Harvey's](https://link.springer.com/article/10.1007/s002140050309) Fortran code that greatly simplifies the setup and includes small improvements gathered from other MECP variants.


## Features

- Single-file Python script. No installation required.
- No need to prepare multiple files that are actually parts of a Gaussian input. Use a normal file with the different values marked and it will just work.
- You don't have to modify the source code anymore. All the setup steps are performed automatically. For example, the number of atoms is inferred from the geometry input.
- Geometry can be specified with atomic symbols or numbers. The script will take care of the rest.
- Convergence criteria can be easily modified.
- No hardcoded values: use another Gaussian version, Fortran compiler, flags...
- Self-contained Python executable with no 3rd party dependencies: you can move it around!
- It still uses the original MECP code behind the scenes, so you can trust the results.
- Optimization trajectory is written for every step.
- Several energy parsers included: `dft`, `pm2`, `ci`, `td`. You can include a new one in a separate Python file if needed.
- Unit-tests (needs `pytest`, `numpy`).

## Why?!

The original MECP code requires a rather intricate setup process. For example, the user is expected to manually modify the code and recompile the Fortran binaries for each calculation. My opinion is that it should be easier and that's why I automated the manual steps in a script called `easymecp`.

## Installation

If you just want to run a few calculations, simply grab the latest `easymecp.py` script from [here](https://github.com/jaimergp/easymecp/blob/master/easymecp/easymecp.py).

If you'd like a permanent installation, use `pip` like this:

```
pip install https://github.com/jaimergp/easymecp/archive/master.zip
```

### Requirements

- Python 2.7/3.4+
- `gfortran` or equivalent compiler
- `pytest` (only for running tests)


## Usage

The `easymecp` package includes a self-contained, pure-Python module called `easymecp` that can be called directly with `python easymecp.py`. If the package is installed with `pip`, the executable `easymecp` will also be available. Either way, the `-h` flag will show the usage:

```
# Direct execution
python easymecp.py -h
# Only available if installed with pip
easymecp -h
```

The MECP program relies on two different input files, one for each multiplicity, segmented in chunks. However, they rarely differ in more than that the spin state, title and `%chk`. EasyMECP supports reading a single Gaussian input file with:

```
python easymecp.py -f system.gjf
```

This file is a normal Gaussian input file, except that the divergent values (multiplicity, chk and title, usually) must be surrounded by curly braces and comma-separated (no spaces).

```
%mem=60GB
%nproc=16
%chk={A,B}.chk
#n PBE1PBE/genecp force guess(read)

{Singlet,Triplet} State

1 {1,3}
   C    0.00000000    1.29390795    0.00000000
   C    1.26894294    0.69491769    0.00000000
   C    1.25319073   -0.68704191    0.00000000
   C    0.00000000   -1.37652121    0.00000000
   C   -1.25319073   -0.68704191    0.00000000
   C   -1.26894294    0.69491769    0.00000000
   H    2.18616903    1.27569862    0.00000000
   H    2.18025862   -1.25246238    0.00000000
   H    0.00000000   -2.46323577    0.00000000
   H   -2.18025862   -1.25246238    0.00000000
   H   -2.18616903    1.27569862    0.00000000
```

Comment lines (normally at the top, but can be anywhere in the file) will be scanned for possible config values if they match this syntax (semicolon can be omitted, spaces are not required around equal sign):

```
! easymecp: max_steps=50
! easymecp TGMax = 7.d-4
```

Please note that if you need an `ExtraOverlays` section before the title, this method would not work. Use the original MECP workflow (explained below) in that case. Individual files for the sections will be automatically generated, so you can use the compatibility mode for easier restarts.

### Output

All the output files will be written to the working directory and stored in the `JOBS` folder, except the one corresponding to the last step. That folder will also contain a `trajectory.xyz` step with all the intermediate geometries. As in the original MECP, `ReportFile` will list more information about the optimization process

### Examples

Run a newstyle MECP calculation with a slightly modified Gaussian input file ([example file](tests/data/C6H5+_singlefile/input.gjf)):

```
python easymecp.py -f input.gjf
```

The file should look like this:

```
! easymecp: max_steps=100
%mem=6GB
%nproc=4
%chk={singlet,triplet}.chk
#n B3LYP/6-31G** force guess(read)

{First,Second} State

1 {1,3}
   6    0.00000000    1.29390795    0.00000000
   6    1.26894294    0.69491769    0.00000000
   6    1.25319073   -0.68704191    0.00000000
   6    0.00000000   -1.37652121    0.00000000
   6   -1.25319073   -0.68704191    0.00000000
   6   -1.26894294    0.69491769    0.00000000
   1    2.18616903    1.27569862    0.00000000
   1    2.18025862   -1.25246238    0.00000000
   1    0.00000000   -2.46323577    0.00000000
   1   -2.18025862   -1.25246238    0.00000000
   1   -2.18616903    1.27569862    0.00000000

```

## Backward compatibility with old MECP workflow

For backward compatibility, a separate mode is provided that mimics the original MECP approach. This mode needs the Gaussian input file separated into different individual files:

- `Input_Header_A`: File containing the top part of system configuration with multiplicity A. Name can be changed with `--a_header` key.
- `Input_Header_B`: File containing the top part of system configuration with multiplicity B. Name can be changed with `--b_header` key.
- `geom`: File containing the starting system geometry (element symbols will be converted in atomic numbers automatically). Name can be changed with `--geom` key.
- `footer` (optional): File containing the bottom part of system configuration (basis sets, etc). Name can be changed with `--footer` flag.

If you use the default filenames and values, you can simply type and the script will run:

```
python easymecp.py
```

Custom configuration values can be specified with command-line flags (for example, `max_steps`).

```
python easymecp.py --geom initial_geometry --FC ifort
```

__TIP__: Use `--gaussian_exe` key to specify the Gaussian executable (version) to use: `g09` or `g16`. Others might work as well. If not specified, `easymecp` will use the default value: `g16` if Gaussian 16 is present in `$PATH`; `g09` otherwise. You can also specify absolute paths here, if needed.

Alternatively, all flags can be written inside a special `*.conf` file (better for long commands and reproducibility) and the program started with:

```
python easymecp.py --conf my.conf
```

More details can be found in the command line help message with `python easymecp.py -h` or `easymecp -h`.


### Old-style examples

Run an old-style MECP calculation providing all files with default names in the same directory as `easymecp.py`:

```
python easymecp.py
```

Increase the maximum allowed steps:

```
python easymecp.py --max_steps 100
```

Less tight TDE convergence criteria (you must use Fortran-style double-precision floats)

```
python easymecp.py --TDE 5.d-4
```

Provide a different energy parser (example file in `tests/data/energy.py`; must include a `parse_energy` top-level function)

```
python easymecp.py --energy_parser energy.py
```

### Restarting jobs

If the calculation does not converge before reaching `max_steps`, you might want to extend it. It's simple: take the latest geometry you like (last one is always in the file `geom`, but you can create another if you want), and either:

A. Manually modify your ``input.gjf`` file... OR

B. Use the old-style options and relaunch with `easymecp --geom <your_geometry>`. Since the new-style will have created the individuals files anyway, you can combine both seamlessly.

No files are lost with the restarts, if the `JOBS` folder already exists, new ones will be automatically named as `JOBS1`, `JOBS2`, etc.

# Cite this work

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4293421.svg)](https://doi.org/10.5281/zenodo.4293421)

```
@software{easymecp,
  author       = {Jaime Rodríguez-Guerra, and Ignacio Funes-Ardoiz, and Feliu Maseras},
  title        = {EasyMECP},
  month        = nov,
  year         = 2018,
  publisher    = {Zenodo},
  version      = {v0.3.2},
  doi          = {10.5281/zenodo.4293421},
  url          = {https://doi.org/10.5281/zenodo.4293421}
}
```


# More MECP implementations for Gaussian

* [MECPro](http://jur.byu.edu/?p=22227): A MIT-licensed, Python-only fork of MECP, developed by J. Snyder and D. H. Ess at BYU.
* [MECPy](http://www2.chemia.uj.edu.pl/~mradon/mecpy/): Dr. Mariusz Radoń developes a separate, Python-only MECP implementation.
* [sobMECP](http://sobereva.com/286). Slight modifications to the original MECP Fortran code and shell scripts so it's a bit automated.
