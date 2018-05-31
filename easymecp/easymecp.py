#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
easyMECP
--------

easymecp v<VERSION>
By: Jaime Rodríguez-Guerra <@jaimergp>, Ignacio Funes

Simplified Python wrapper around the original MECP Fortran code
from J. Harvey (2003).

It still uses the same Fortran code*, but saves the hassle of manual
edition of files** and shell files.

* gfortran or equivalent is still required behind the scenes
** Number of atoms is now inferred from geometry file automatically

Recommended workflow (-f)
.........................

The MECP program expects two different input files, one for each
multiplicity. However, they rarely differ in more than that. EasyMECP
supports reading a single Gaussian input file with

    python easymecp.py -f system.gjf

This file is a normal Gaussian input file, except that the divergent
values (multiplicity, chk and title, usually) must be surrounded by
curly braces and comma-separated (no spaces).

    %mem=60GB
    %nproc=16
    %chk={A,B}.chk
    #n PBE1PBE/genecp force guess(read)

    {Singlet,Triplet} State

    1 {1,3}

Comment lines (normally at the top, but can be anywhere in the file)
will be scanned for possible config values if they match this syntax
(semicolon can be omitted, spaces are not required around equal sign):

    ! easymecp: max_steps=50
    ! easymecp TGMax = 7.d-4

Please note that if you need ExtraOverlays, this method would not work.
Use the original MECP workflow (explained below) in that case.

Individual files for the sections will be automatically generated, so you
can use the compatibility mode for easier restarts.

Compatibility mode
..................

If you prefer to use the original MECP workflow, you can also specify
individual files:

- `Input_Header_A`: Input lines before geometry (route, %-lines, title, charge,
  and spin for multiplicity A) . You can specify a custom filename with
  --header_a flag.
- `Input_Header_B`: Input lines before geometry (route, %-lines, title, charge,
  and spin for multiplicity B) . You can specify a custom filename with
  --header_b flag.
- `geom`: Initial geometry (choose another name with --geom)
- `footer`: bottom contents (basis sets, modredundant, etc). Choose another name
   with --footer, if desired.


These options (and more) can be specified in the command line with the
appropriate flags (see easymecp -h for the full list). For example

    python easymecp.py --geom initial_geometry

If the command gets too long, you can opt to use a config file that
specifies key-value pairs like this:

    geom: initial_geometry
    footer: filetail

And them have it called like this:

    python easymecp.py --conf my.conf

Command line arguments, input and config files CANNOT be mixed. Input files
have higher precedence over config files, and this one over command-line arguments.

    -f > --conf > anything else

"""

from __future__ import print_function, absolute_import
import sys
if sys.version_info[0:2] < (3, 4) and sys.version_info[0:2] != (2, 7):
    sys.exit('! ERROR: easyMECP requires Python 2.7 or 3.4+')

from contextlib import contextmanager
from datetime import datetime
from distutils.spawn import find_executable
from runpy import run_path
from subprocess import call, SubprocessError
from tempfile import mkdtemp
import argparse
import os
import re
import shlex
import shutil


__version__ = '0.2.0'


class MECPCalculation(object):

    OK = 'OK'
    ERROR = 'ERROR'
    MAX_ITERATIONS_REACHED = 'MAX_ITERATIONS_REACHED'

    def __init__(self, max_steps=50, a_header='Input_Header_A', b_header='Input_Header_B',
                 geom='geom', footer='footer', natom=0, TDE='5.d-5',
                 TDXMax='4.d-3', TDXRMS='2.5d-3', TGMax='7.d-4', TGRMS='5.d-4',
                 energy_parser='dft', FC=os.environ.get('FC', 'gfortran'),
                 FFLAGS=os.environ.get('FFLAGS', '-O -ffixed-line-length-none'),
                 gaussian_exe=('g16' if find_executable('g16') else 'g09'), **kwargs):
        self.a_header = _extant_file(a_header, name='a_header')
        self.b_header = _extant_file(b_header, name='b_header')
        self.geom = _extant_file(geom, name='geom')
        self.footer = _extant_file(footer, name='footer', allow_errors=True)
        self.max_steps = max_steps
        self.gaussian_exe = gaussian_exe
        self.TDE = TDE
        self.TDXMax = TDXMax
        self.TDXRMS = TDXRMS
        self.TGMax = TGMax
        self.TGRMS = TGRMS
        self.FC = FC
        self.FFLAGS = FFLAGS
        self.natom = natom
        self.converged_at = None
        if energy_parser.endswith('.py') and os.path.isfile(energy_parser):
            try:
                self._parse_energy = run_path(energy_parser)['parse_energy']
            except KeyError:
                raise ValueError('External energy parsers must define a `parse_energy` function')
        elif callable(energy_parser):
            self._parse_energy = energy_parser
        else:
            try:
                self._parse_energy = globals()['_parse_energy_' + energy_parser]
            except KeyError:
                raise ValueError('energy_parser `{}` must be one of <{}>'.format(
                                energy_parser, ', '.join(AVAILABLE_ENERGY_PARSERS)))
        if not natom:
            with open(geom) as f:
                for line in f:
                    if line.startswith('!'):
                        continue
                    elif len(line.split()) >= 4:
                        self.natom +=1

    @classmethod
    def from_conf(cls, path):
        d = _get_defaults()
        with open(path) as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#') or ':' not in line:
                    continue
                fields = line.split(':', 1)
                if len(fields) < 2 or not fields[1].strip():
                    print('Malformed line #{}: {}'.format(i, line))
                    sys.exit()
                key, value = fields[0], fields[1].strip()
                if key not in DEFAULTS:
                    print('! Skipping key `{}` (not recognized)'.format(key))
                    continue
                if key in ('natom', 'max_steps'):
                    value = int(value)
                elif key in ('TDE', 'TDXMax', 'TDXRMS', 'TGMax', 'TGRMS'):
                    value = _valid_fortran_double(value, key)
                elif key in ('a_header', 'b_header', 'footer', 'geom'):
                    if not os.path.isfile(value):
                        print('! `{}` file with path `{}` not available!'.format(key, value))
                        sys.exit()
                d[key] = value
        return cls(**d)

    @classmethod
    def from_gaussian_input_file(cls, path):
        def _process_header_line(line):
            match = re.search(r'{(.*),(.*)}', line)
            if match:
                return (line.replace(match.group(0), match.group(1)),
                        line.replace(match.group(0), match.group(2)))
            return line, line

        d = _get_defaults()
        section = 0
        header_a, header_b, geom, footer = [], [], None, []
        # Parse Gaussian input files into header(s), geometry and footer
        with open(path) as f:
            for line in f:
                # Detect sections
                if not line.strip():
                    section += 1

                # Assign lines to sections
                # <HEADER>
                if section <= 1:
                    line_a, line_b = _process_header_line(line)
                    header_a.append(line_a)
                    header_b.append(line_b)
                elif section == 2:
                    if not line.strip():
                        header_a.append(line)
                        header_b.append(line)
                    elif geom is None:
                        line_a, line_b = _process_header_line(line)
                        header_a.append(line_a)
                        header_b.append(line_b)
                        geom = []
                # </HEADER>
                # <GEOM>
                    else:
                        # Everything above is part of the header
                        geom.append(line)
                # </GEOM>
                # <FOOTER>
                elif section >= 3:
                    footer.append(line)
                #</FOOTER>

                # Detect special comments
                if line.startswith('!'):
                    match = re.search(r'^! easymecp:?\s+(\S+)\s*=(\S+)', line)
                    if match:
                        key = match.group(1)
                        value = match.group(2)
                        if key == 'max_steps':
                            d[key] = int(value)
                        elif key in ('TDE', 'TDXMax', 'TDXRMS', 'TGMax', 'TGRMS'):
                            d[key] = _valid_fortran_double(value, key)

        # Write temporary files
        d['a_header'] = '_header_a'
        d['b_header'] = '_header_b'
        d['geom'] = '_initial_geom'
        d['footer'] = '_footer'
        for filepath, lines in (('_header_a', header_a), ('_header_b', header_b),
                                ('_initial_geom', geom), ('_footer', footer)):
            if lines:
                with open(filepath, 'w') as f:
                    f.writelines(lines)
        with open(os.path.splitext(os.path.basename(path))[0] + '.conf', 'w') as f:
            f.write('\n'.join('{}: {}'.format(k, v) for (k,v) in d.items()))

        return cls(**d)

    def run(self):
        # Recompile Fortran code
        print('Running easyMECP v{}...'.format(__version__))
        print('Preparing workspace...')
        self.prepare_workspace()
        print('Compiling MECP for {} atoms...'.format(self.natom))
        mecp_exe = self.compile_fortran()
        geom = self.geom

        for i in range(self.max_steps):
            # First, compute and parse Gaussian Jobs
            print('Running step #{}...'.format(i))
            self.add_trajectory_step(geom, step=i)

            print('  Launching Gaussian job for file A...')
            input_a = self.prepare_gaussian(self.a_header, geom, self.footer, label='A', step=i)
            energy_a, gradients_a = self.run_gaussian(input_a)
            print('  Launching Gaussian job for file B...')
            input_b = self.prepare_gaussian(self.b_header, geom, self.footer, label='B', step=i)
            energy_b, gradients_b = self.run_gaussian(input_b)

            # Second, run MECP
            if not self.prepare_ab_initio(energy_a, energy_b, gradients_a, gradients_b):
                return self.ERROR
            print('  Launching MECP...')
            try:
                retcode = call([mecp_exe], stdout=sys.stdout, stderr=sys.stderr)
                if retcode:
                    raise SubprocessError('MECP returned code {}'.format(retcode))
            except Exception as e:
                print('  ! Error during MECP execution:', e.__class__.__name__, '->', e)
                self.report('ERROR')
            if os.path.isfile('AddtoReportFile'):
                with open('AddtoReportFile') as f:
                    self.report(f.read())
                os.remove('AddtoReportFile')

            # Third, check results in this iteration
            with open('ReportFile') as f:
                for line in f:
                    if 'CONVERGED' in line:
                        print('  MECP optimization has converged at Step', i)
                        self.converged_at = i
                        return self.OK
                    if 'ERROR' in line:
                        print('  An error ocurred!')
                        return self.ERROR
            # Keep going! Next iteration will use geom autogenerated by MECP.x
            geom = 'geom'
            os.rename('Job{}_A.gjf'.format(i), 'JOBS/Job{}_A.gjf'.format(i))
            os.rename('Job{}_A.log'.format(i), 'JOBS/Job{}_A.log'.format(i))
            os.rename('Job{}_B.gjf'.format(i), 'JOBS/Job{}_B.gjf'.format(i))
            os.rename('Job{}_B.log'.format(i), 'JOBS/Job{}_B.log'.format(i))
            print()

        return self.MAX_ITERATIONS_REACHED

    def compile_fortran(self):
        with temporary_directory(enter=False) as tmp:
            # Patch source code
            code = MECP_FORTRAN.format(NUMATOM=self.natom, TDE=self.TDE, TDXMax=self.TDXMax,
                                       TDXRMS=self.TDXRMS, TGMax=self.TGMax, TGRMS=self.TGRMS)
            with open(os.path.join(tmp, 'MECP.f'), 'w') as f:
                f.write(code)
            with open('_fortran_compilation', 'w') as out:
                call([self.FC] + shlex.split(self.FFLAGS) + ['MECP.f', '-o', 'MECP.x'], cwd=tmp, stdout=out)
            shutil.copyfile(os.path.join(tmp, 'MECP.x'), 'MECP.x')
        os.chmod('MECP.x', os.stat('MECP.x').st_mode | 0o111)  # make executable
        return './MECP.x'

    def element_symbol_to_number(self, fh, drop_blank=True):
        elements = ELEMENTS
        lines = []
        for line in fh:
            if drop_blank and not line.strip():
                continue
            fields = line.split()
            if len(fields) > 3:
                if not fields[0].isdigit():
                    number = elements.get(fields[0].title(), -1)
                    line = line.replace(fields[0], str(number))
            lines.append(line)
        return ''.join(lines)

    def element_number_to_symbol(self, fh, drop_blank=False):
        elements = REVERSE_ELEMENTS
        lines = []
        for line in fh:
            if drop_blank and not line.strip():
                continue
            fields = line.split()
            if len(fields) > 3:
                if fields[0].isdigit():
                    symbol = elements.get(int(fields[0]), 'LP')
                    line = line.replace(fields[0], symbol, 1)
            lines.append(line)
        return ''.join(lines)

    def prepare_workspace(self):
        i = 1
        while os.path.exists('JOBS'):
            try:
                os.rename('JOBS', 'JOBS{}'.format(i))
            except OSError:
                i += 1
        os.makedirs('JOBS')
        with open(self.geom) as f:
            geometry = self.element_symbol_to_number(f)
        with open('ProgFile', 'w') as f:
            f.seek(0)
            f.write(PROGFILE.format(natom=self.natom, geometry=geometry))
            f.truncate()
        with open('ReportFile', 'w') as f:
            f.seek(0)
            f.write('{}\n'.format(datetime.now()))
            f.truncate()

    def prepare_gaussian(self, header, geom, footer, label='A', step=0):
        name = 'Job{}_{}.gjf'.format(step, label)
        with open(name, 'w') as f:
            with open(header) as a:
                contents = a.read()
                chk = re.search(r'^%chk=(.*)$', contents, flags=re.IGNORECASE|re.MULTILINE)
                chk_exists = chk and chk.group(1) and os.path.isfile(chk.group(1).strip())
                if not step and not chk_exists:
                    contents = contents.replace('guess(read)', '')
                f.write(contents.rstrip())
            f.write('\n')
            with open(geom) as b:
                contents = self.element_symbol_to_number(b, drop_blank=True)
                f.write(contents.rstrip())
            f.write('\n\n')
            try:
                with open(footer) as c:
                    f.write(c.read().lstrip())
            except IOError:
               pass
            f.write('\n\n')  # Gaussian is picky about file endings...
        return name

    def run_gaussian(self, inputfile):
        try:
            call([self.gaussian_exe, inputfile], stdout=sys.stdout, stderr=sys.stderr)
        except Exception as e:
            print('  ! Could not run Gaussian job', inputfile)
            print('  !', e.__class__.__name__, '->', e)
            self.report('ERROR')
            return None, None
        energy = None
        gradients = []
        with open(os.path.splitext(inputfile)[0] + '.log') as f:
            for line in f:
                fields = line.split()
                if 'Convergence failure' in line:
                    print('  ! There has been a convergence problem with', inputfile)
                    self.report('ERROR')
                    return None, None
                # gradients
                gradients = self._parse_gradients(f, line, fields, default=gradients)
                # energies
                energy = self._parse_energy(f, line, fields, default=energy)
        return energy, gradients

    def _parse_gradients(self, f, line, fields, default=None):
        if len(fields) > 2 and fields[0] == 'Center' and fields[1] == 'Atomic' and fields[2] == 'Forces':
            gradients = []
            next(f), next(f)  # skip two lines
            for _ in range(self.natom):
                line = next(f)
                gradients.append(line.split()[1:5])
            return gradients
        return default

    def prepare_ab_initio(self, energy_a, energy_b, gradients_a, gradients_b):
        if not all(x is not None for x in (energy_a, energy_b, gradients_a, gradients_b)):
            print('  ! Some energies or gradients could not be obtained!')
            return False
        with open('ab_initio', 'w') as f:
            print('Energy of the First State', file=f)
            print(energy_a, file=f)
            print('Gradient of the First State' , file=f)
            print('\n'.join(['{}  {}   {}   {}'.format(*l) for l in gradients_a]), file=f)
            print('Energy of the Second State' , file=f)
            print(energy_b, file=f)
            print('Gradient of the Second State' , file=f)
            print('\n'.join(['{}  {}   {}   {}'.format(*l) for l in gradients_b]), file=f)
        return True

    def report(self, msg):
        with open('ReportFile', 'a') as r:
            print(msg, file=r)

    def add_trajectory_step(self, geometry='geom', step=0):
        with open('JOBS/trajectory.xyz', 'a') as f:
            print(self.natom, file=f)
            print('Step', step, file=f)
            with open(geometry) as g:
                print(self.element_number_to_symbol(g, drop_blank=True), file=f)


###
# Energy parsers
###
def _parse_energy_dft(f, line, fields, default=None):
    if len(fields) > 4 and fields[0] == 'SCF' and fields[1] == 'Done:':
        return float(fields[4])
    return default


def _parse_energy_mp2(f, line, fields, default=None):
    if len(fields) > 5 and fields[0] == 'E2' and fields[3] == 'EUMP2':
        return float(fields[5])
    return default


def _parse_energy_cis(f, line, fields, default=None):
    if len(fields) > 4 and fields[2] == 'E(CIS)':
        return float(fields[4])
    return default


def _parse_energy_td(f, line, fields, default=None):
    if len(fields) > 4 and fields[2] == 'E(TD-HF/TD-KS)':
        return float(fields[4])
    return default


###
# Validators
###
def _valid_fortran_double(value, key=None):
    try:
        float(value.replace('d', 'e'))
    except ValueError:
        raise ValueError('Value {}{} must be a valid Fortran double ([base].d[exp], like 5.d-4)'.format(
                         value, 'for {}'.format(key) if key is not None else ''))
    else:
        return value


def _extant_file(path, name=None, allow_errors=False):
    """ Verify file exists or report error """
    if os.path.isfile(path):
        return path
    label = '`{}` with path '.format(name) if name is not None else ''
    skipping = ' Skipping...' if allow_errors else ''
    msg = 'File {label}{path} is not available!{skipping}'.format(
          label=label, path=path, skipping=skipping)
    if allow_errors:
        print(msg)
        return path
    else:
        raise ValueError(msg)


###
# Helpers
###
def _get_defaults():
    _defaults = MECPCalculation.__init__.__defaults__
    _ndef = len(_defaults)
    _narg = MECPCalculation.__init__.__code__.co_argcount
    _args = MECPCalculation.__init__.__code__.co_varnames[(_narg-_ndef):_narg]
    return dict(zip(_args, _defaults))


@contextmanager
def temporary_directory(enter=True, **kwargs):
    """Create and enter a temporary directory; used as context manager."""
    temp_dir = mkdtemp(**kwargs)
    if enter:
        cwd = os.getcwd()
        os.chdir(temp_dir)
    yield temp_dir
    if enter:
        os.chdir(cwd)
    shutil.rmtree(temp_dir)


###
# App
###
def _parse_cli():
    description = __doc__.replace('<VERSION>', __version__)
    p = argparse.ArgumentParser(prog='easymecp', description=description,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-V', '--version', action='version', version='%(prog)s v' + __version__)
    p.add_argument('-f', '--inputfile', metavar='INPUTFILE', help='Initialize from Gaussian input file. '
    'Divergent options can be specified with curly braces at any time: {A,B}. Additional flags '
    'must be specified in comment lines (read above).')
    p.add_argument('--conf', metavar='CONFFILE', help='Initialize from configuration file. '
                   'Each value must be provided in its own line, with syntax <key>: <value>.')
    defaults = _get_defaults()
    for k, v in sorted(defaults.items()):
        p.add_argument('--'+k, default=v, metavar='VALUE', help='{} (default={!r})'.format(USAGE[k], v))
    args = p.parse_args()
    return args


def main():
    args = _parse_cli()
    try:
        if args.inputfile:
            print('Gaussian input file', args.inputfile, 'has been specified. '
                  'Ignoring any other arguments!')
            calc = MECPCalculation.from_gaussian_input_file(args.inputfile)
        elif args.conf:
            print('Configuration file', args.conf, 'has been specified. '
                  'Ignoring any other arguments!')
            calc = MECPCalculation.from_conf(args.conf)
        else:
            calc = MECPCalculation(**args.__dict__)
    except ValueError as e:
        print('! ERROR:', e)
        print('         Run easymecp -h (or python easymecp.py -h) for help.')
        sys.exit()
    result = calc.run()
    if result == MECPCalculation.OK:
        print('Success! Check ReportFile for results.')
    elif result == MECPCalculation.ERROR:
        print('Something failed... Check Gaussian outputs and/or ReportFile.')
    elif result == MECPCalculation.MAX_ITERATIONS_REACHED:
        print('Calculation did not converge after', calc.max_steps, 'steps...')


###
# Constants
###
DEFAULTS = _get_defaults()
AVAILABLE_ENERGY_PARSERS = set([key[14:] for key in globals().copy()
                                if key.startswith('_parse_energy_')])
PROGFILE = """
Title
Number of Atoms
{natom}
Number of Steps already Run
0
Is this a complete ProgFile ?
0
Geometry:
{geometry}
"""
ELEMENTS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
    'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17,
    'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
    'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33,
    'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41,
    'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
    'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
    'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
    'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
    'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81,
    'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89,
    'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
    'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104,
    'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111,
    'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}
REVERSE_ELEMENTS = dict((v, k) for (k, v) in ELEMENTS.items())
USAGE = {
    'max_steps': 'Number of iterations to perform until stop or convergence',
    'natom': 'If easymecp cannot detect number of atoms automatically, override with this option',
    'a_header': 'File containing the top part of system configuration with multiplicity A',
    'b_header': 'File containing the top part of system configuration with multiplicity B',
    'geom': 'File containing the starting system geometry '
            '(element symbols will be converted in atomic numbers automatically)',
    'footer': 'File containing the bottom part of system configuration',
    'energy_parser': 'Which energy should be parsed: dft, mp2, cis, td. It can also be a '
                     'path to a Python file containing a `parse_energy` function.',
    'gaussian_exe': 'Path to gaussian executable. Compatible versions: g09, g16',
    'TDE': 'Convergence threshold for difference in E. Must be a valid Fortran double!',
    'TDXMax': 'Convergence threshold for max change of X. Must be a valid Fortran double!',
    'TDXRMS': 'Convergence threshold for RMS change of X. Must be a valid Fortran double!',
    'TGMax': 'Convergence threshold for max gradient el. Must be a valid Fortran double!',
    'TGRMS': 'Convergence threshold for rms gradient el. Must be a valid Fortran double!',
    'FC': 'Fortran compiler (can also be set with $FC environment variable)',
    'FFLAGS': 'Fortran compiler flags (can also be set with $FFLAGS environment variable)'
}

MECP_FORTRAN = """
      PROGRAM Optimizer
      implicit none

C     New version, Nov. 2003. All gradients converted to Hartree/Ang.
C     facPp chosen so that inverse Hessian ~ diagonal 0.7 Ang*Ang/Hartree
C     Packaged, Apr 2018. Jaime RGP.

      INTEGER Natom, Nx
      parameter (Natom={NUMATOM})
C         Set Natom to the number of atoms.
      parameter (Nx=3*Natom)

      INTEGER Nstep, FFile, Conv, AtNum(Natom)
      DOUBLE PRECISION Ea_1, Ea_2, Eb_1, Eb_2, IH_1(Nx,Nx), IH_2(Nx,Nx), ParG(Nx), PerpG(Nx)
      DOUBLE PRECISION Ga_1(Nx), Ga_2(Nx), Gb_1(Nx), Gb_2(Nx), X_1(Nx), X_2(Nx), X_3(Nx)
      DOUBLE PRECISION G_1(Nx), G_2(Nx)

      CALL ReadProgFile(Natom,Nx,AtNum,Nstep,FFile,X_1,X_2,IH_1,Ea_1,Eb_1,Ga_1,Gb_1,G_1)
C          Recover data from the previous Steps
      CALL ReadInput(Natom,Nx,Ea_2,Ga_2,Eb_2,Gb_2)
C          Read in the new ab initio Energies and Gradients
      CALL Effective_Gradient(Nx,Ea_2,Eb_2,Ga_2,Gb_2,ParG,PerpG,G_2)
C          Compute the Effective Gradient fac * DeltaE * PerpG + ParG
      CALL UpdateX(Nx,Nstep,FFile,X_1,X_2,X_3,G_1,G_2,IH_1,IH_2)
C          the BFGS step
      CALL TestConvergence(Nx,Natom,Nstep,AtNum,Ea_2,Eb_2,X_2,X_3,ParG,PerpG,G_2,Conv)
C          Checks Delta E, X, Mag. Delta G. Writes Output File
      IF (Conv .ne. 1) THEN
           CALL WriteGeomFile(Natom,Nx,AtNum,X_3)
           CALL WriteProgFile(Natom,Nx,AtNum,Nstep,X_2,X_3,IH_2,Ea_2,Eb_2,Ga_2,Gb_2,G_2)
      END IF

      END

       SUBROUTINE Effective_Gradient(N,Ea,Eb,Ga,Gb,ParG,PerpG,G)
       implicit none

C      Computes the parallel and perpendicular compenents of the Effective Gradient,
C      As well as the effective gradient itself.

       integer n, i
       double precision Ea, Eb, Ga(N), Gb(N), G(N), ParG(N), PerpG(N), npg, pp
       double precision facPP, facP
       parameter (facPP=140.d0,facP=1.d0)
C  These factors are only really important for the first step
C  The "difference gradient" is typically ca. 0.075 Hartree/Bohr.
C      i.e. 0.14 Hartree/Angstrom.
C  Assuming this is constant, this gives a Hessian term for the func (Ea-Eb)**2
C     of ca. 0.01 Hartree**2 / Angstrom**2  (0.14**2 / 2)
C  The Hessian term along "normal" coordinates is empirically about 1.4 Hartree / Angstrom**2
C  Using facPP ~ 140 /Hartree means that along this coordinate, too, the Hessian is about right.

       npg = 0.d0
       pp = 0.d0
       do i = 1, n
           PerpG(i) = Ga(i) - Gb(i)
           npg = npg + PerpG(i)**2
           pp = pp + Ga(i) * PerpG(i)
       end do
       npg = sqrt(npg)
       pp = pp / npg
       do i = 1, n
           ParG(i) = Ga(i) - PerpG(i) / npg * pp
           G(i) = (Ea - Eb) * facPP * PerpG(i) + facP * ParG(i)
       end do
       return

       END

      SUBROUTINE error_handling(e)
      implicit none

      integer e

      IF (e .eq. -1) THEN
          write (*,*) "The ProgFile does not exist for some reason."
          write (*,*) "Please create it then try again."
      END IF

      IF (e .eq. -2) THEN
          write (*,*) "Incomplete ProgFile - did you set Nstep right ?"
      END IF

      IF (e .eq. -3) THEN
          write (*,*) "Problem Reading Gaussian Output"
          OPEN(UNIT=302,FILE="AddtoReportFile")
          write (302,*) "ERROR IN GAUSSIAN STEP"
          CLOSE(302)
      END IF

      IF (e .eq. -4) THEN
           write (*,*) "Wrong Number of Atoms - Recompile !!"
      END IF

      IF (e .eq. -5) THEN
           write (*,*) "Problem with the ab initio Job."
      END IF

      STOP

      END

       SUBROUTINE Initialize(N,H,Ea,Eb,Ga,Gb)
       implicit none

       INTEGER i, j, n
       DOUBLE PRECISION H(N,N), Ea, Eb, Ga(N), Gb(N)

       Ea = 0.d0
       Eb = 0.d0
       do i = 1, n
           Ga(i) = 0.d0
           Gb(i) = 0.d0
           H(i,i) = .7d0
C Inverse Hessian values of .7 correspond to Hessians of 1.4 Hartree/Angstrom**2 - about right
           do j = i + 1, n
               H(i,j) = 0.d0
               H(j,i) = 0.d0
           end do
       end do
       return

       END

       SUBROUTINE ReadInput(natom,nx,Ea,Ga,Eb,Gb)
       implicit none

       INTEGER i, j, k, kk, natom, nx, error
       DOUBLE PRECISION Ea, Eb, Ga(Nx), Gb(Nx)
       CHARACTER*60 dummy
       DOUBLE PRECISION tmpg, bohr
       PARAMETER (bohr=0.529177d0)

       OPEN(UNIT=8,FILE="ab_initio")

       read (8,*) dummy
       read (UNIT=8,FMT=*,IOSTAT=error) Ea
       IF (error .ne. 0) THEN
            error = -5
            CALL error_handling(error)
       END IF
       read (8,*) dummy
       DO I = 1, natom
           k = 3 * (i - 1) + 1
           READ (UNIT=8,FMT=*,IOSTAT=error) kk, (Ga(J),J=k,k+2)
           IF (error .ne. 0) THEN
                error = -5
                CALL error_handling(error)
           END IF
       END DO
       read (8,*) dummy
C The gradient output by most programs is in fact the FORCE - convert here to gradient.
C Also, convert from Hartree/Bohr to Hartree/Ang
       DO i = 1, nx
           Ga(i) = -Ga(i) / bohr
       end do

       read (UNIT=8,FMT=*,IOSTAT=error) Eb
       read (8,*) dummy
       IF (error .ne. 0) THEN
            error = -5
            CALL error_handling(error)
       END IF
       DO I = 1, natom
           k = 3 * (i - 1) + 1
           READ (UNIT=8,FMT=*,IOSTAT=error) kk, (Gb(J),J=k,k+2)
           IF (error .ne. 0) THEN
                error = -5
                CALL error_handling(error)
           END IF
       END DO
       DO i = 1, nx
           Gb(i) = -Gb(i) / bohr
       end do

       CLOSE(8)
       return

       END
       SUBROUTINE ReadProgFile(natom,nx,AtNum,Nstep,FFile,X_1,X_2,HI,Ea,Eb,Ga,Gb,G)
       implicit none

       INTEGER natom, nx, AtNum(Natom), Nstep, i, j, k, error, FFile
       DOUBLE PRECISION X_1(Nx), X_2(Nx), HI(Nx,Nx), Ea, Eb, Ga(Nx), Gb(Nx), G(Nx)

       LOGICAL PGFok
       CHARACTER*80 dummy

       INQUIRE(FILE="ProgFile",EXIST=PGFok)
       IF (.not. PGFok) THEN
           error = -1
           CALL error_handling(error)
       END IF

       OPEN(UNIT=8,FILE="ProgFile")
       READ(8,*) dummy
       READ(8,*) dummy
       READ(8,*) i
       IF (i .ne. natom) then
           error = -4
           CALL error_handling(error)
       END IF
       READ(8,*) dummy
       READ(8,*) Nstep
       READ(8,*) dummy
       READ(8,*) FFile

       READ(8,*) dummy
       DO I = 1, natom
           k = 3 * (i - 1) + 1
           READ(8,*) AtNum(I), (X_2(J),J=k,k+2)
       END DO

       IF ((Nstep .eq. 0) .AND. (FFile .eq. 0)) THEN
C First step - the ProgFile is incomplete. To avoid problems,
C there are 'existence' tests in the next input lines so
C as to crash the program if Nstep > 0 yet the file is incomplete.
           CALL Initialize(Nx,HI,Ea,Eb,Ga,Gb)
           CLOSE(8)
           RETURN
       END IF

       READ(UNIT=8,FMT=*,IOSTAT=error) dummy
       IF (error .ne. 0) THEN
           error = -2
           CALL error_handling(error)
       END IF
       DO I = 1, Natom
           k = 3 *(I-1) + 1
           READ(UNIT=8,FMT=*,IOSTAT=error) (X_1(J),J=k,k+2)
C Geometries are in Angstrom - fine!
           IF (error .ne. 0) THEN
               error = -2
               CALL error_handling(error)
           END IF
       END DO
       READ(8,*) dummy
       READ(8,*) Ea
       READ(8,*) Eb
       READ(8,*) dummy
       DO I = 1, nx
           READ(8,*) Ga(I)
       END DO
       READ(8,*) dummy
       DO I = 1, nx
           READ(8,*) Gb(I)
       END DO
       READ(8,*) dummy
       DO I = 1, nx
           READ(8,*) G(I)
       END DO
       READ(8,*) dummy
       DO I = 1, nx
           DO J = 1, nx
               READ(8,*) HI(I,J)
           END DO
       END DO
       CLOSE(8)
       return

       END
      SUBROUTINE TestConvergence(N,Natom,Nstep,AtNum,Ea,Eb,X_2,X_3,ParG,PerpG,G,Conv)
      implicit none

C     Checks convergence, and updates report file
C     There are 5 criteria for testing convergence
C     They are the same as in Gaussian (Except the last):
C     Av.DeltaX, Max.DeltaX, Av.Grad., Max.Grad., DeltaE

      INTEGER N, Natom, AtNum(Natom), Nstep, Conv, i, j, k
      DOUBLE PRECISION Ea, Eb, X_2(N), ParG(N), PerpG(N), X_3(N), G(N)

      CHARACTER*3 flags(5)
      LOGICAL PConv(5)
      DOUBLE PRECISION DeltaX(N), DE, DXMax, DXRMS, GMax, GRMS, PpGRMS, PGRMS, TDE
      DOUBLE PRECISION TDXMax, TDxRMS, TGMax, TGRMS
      PARAMETER (TDE={TDE},TDXMax={TDXMax},TDXRMS={TDXRMS},TGMax={TGMax},TGRMS={TGRMS})

      OPEN(UNIT=8,FILE="AddtoReportFile")
      IF (Nstep .eq. 0) THEN
          write (8,*) "      Geometry Optimization of an MECP"
          write (8,*) "      Program: J. N. Harvey, March 1999"
          write (8,*) "        version 2, November 2003"
          write (8,*) "      easyMECP: J. RG. Pedregal, May 2018"
          write (8,*)
          write (8,'(A)') "Initial Geometry:"
          DO I = 1, Natom
              k = 3 * (I-1) + 1
              write (8,'(I3,3F15.7)') AtNum(I), (X_2(j),j=k,k+2)
          END DO
          write (8,*)
      END IF
      DE = ABS(Ea - Eb)
      DXMax = 0.d0
      DXRMS = 0.d0
      GMax = 0.d0
      GRMS = 0.d0
      PpGRMS = 0.d0
      PGRMS = 0.d0
      DO I = 1, n
          DeltaX(i) = X_3(i) - X_2(i)
          IF (ABS(DeltaX(i)) .gt. DXMax) DXMax = ABS(DeltaX(i))
          DXRMS = DXRMS + DeltaX(i)**2
          IF (ABS(G(i)) .gt. Gmax) Gmax = ABS(G(i))
          GRMS = GRMS + G(i)**2
          PpGRMS = PpGRMS + PerpG(i)**2
          PGRMS = PGRMS + ParG(i)**2
      END DO
      DXRMS = SQRT(DXRMS / N)
      GRMS = SQRT(GRMS / N)
      PpGRMS= SQRT(PpGRMS / N)
      PGRMS = SQRT(PGRMS / N)

      Conv = 0
      do i = 1, 5
          flags(i) = " NO"
          PConv(i) = .false.
      end do

      IF (GMax .lt. TGMax) THEN
          PConv(1) = .true.
          flags(1) = "YES"
      END IF
      IF (GRMS .lt. TGRMS) THEN
          PConv(2) = .true.
          flags(2) = "YES"
      END IF
      IF (DXMax .lt. TDXMax) THEN
          PConv(3) = .true.
          flags(3) = "YES"
      END IF
      IF (DXRMS .lt. TDXRMS) THEN
          PConv(4) = .true.
          flags(4) = "YES"
      END IF
      IF (DE .lt. TDE) THEN
          PConv(5) = .true.
          flags(5) = "YES"
      END IF
      IF (PConv(1) .and. PConv(2) .and. PConv(3) .and. PConv(4) .and. PConv(5)) THEN
          Conv = 1
      ELSE
          Nstep = Nstep + 1
      END IF

      write (8,'(A,F18.10)') "Energy of First State:  ",Ea
      write (8,'(A,F18.10)') "Energy of Second State: ",Eb
      write (8,*)
      write (8,'(A)') "Convergence Check (Actual Value, then Threshold, then Status):"
      write (8,'(A,F11.6,A,F8.6,A,A)') "Max Gradient El.:", GMax," (",TGMax,")  ",flags(1)
      write (8,'(A,F11.6,A,F8.6,A,A)') "RMS Gradient El.:", GRMS," (",TGRMS,")  ",flags(2)
      write (8,'(A,F11.6,A,F8.6,A,A)') "Max Change of X: ",DXMax," (",TDXMax,")  ",flags(3)
      write (8,'(A,F11.6,A,F8.6,A,A)') "RMS Change of X: ",DXRMS," (",TDXRMS,")  ",flags(4)
      write (8,'(A,F11.6,A,F8.6,A,A)') "Difference in E: ",DE," (",TDE,")  ",flags(5)
      write (8,*)
      write (8,'(A)') "Overall Effective Gradient:"
      DO I = 1, Natom
          k = 3 * (I - 1) + 1
          write (8,'(I3,3F16.8)') I, (G(J),J=k,k+2)
      END DO
      write (8,*)
      write (8,'(A,F11.6,A)') "Difference Gradient: (RMS * DE:",PpGRMS,")"
      DO I = 1, Natom
          k = 3 * (I - 1) + 1
          write (8,'(I3,3F16.8)') I, (PerpG(J),J=k,k+2)
      END DO
      write (8,*)
      write (8,'(A,F11.6,A)') "Parallel Gradient: (RMS:",PGRMS,")"
      DO I = 1, Natom
          k = 3 * (I - 1) + 1
          write (8,'(I3,3F16.8)') I, (ParG(J),J=k,k+2)
      END DO
      write (8,*)

      IF (Conv .eq. 1) THEN
          write (8,'(A)') "The MECP Optimization has CONVERGED at that geometry !!!"
          write (8,'(A)') "Goodbye and fly with us again..."
      ELSE
          write (8,'(A,I3)') "Geometry at Step",Nstep
          DO I = 1, Natom
              k = 3 * (I - 1) + 1
              write (8,'(I3,3F15.7)') AtNum(I), (X_3(J),J=k,k+2)
          END DO
          write (8,*)
      END IF

      CLOSE(8)
      return

      END

       SUBROUTINE UpdateX(N,Nstep,FFile,X_1,X_2,X_3,G_1,G_2,HI_1,HI_2)
       implicit none

       ! Specially Adapted BFGS routine from Numerical Recipes

       integer i, j, k, n, Nstep, FFile
       double precision X_1(N), X_2(N), G_1(N), G_2(N), HI_1(N,N), X_3(N), HI_2(N,N)
C change fmc
       double precision stpmax, DelG(N), HDelG(N), ChgeX(N), DelX(N), w(N),
     1 fac,fad, fae, sumdg, sumdx, stpl, lgstst, STPMX
       parameter (STPMX = 0.1d0)

       stpmax = STPMX * N
       IF ((Nstep .eq. 0) .and. (FFile .eq. 0)) THEN
           DO i = 1, n
               ChgeX(i) = -.7d0 * G_2(i)
               DO j = 1, n
                   HI_2(i,j) = HI_1(i,j)
               end do
           end do
       ELSE
           DO i = 1, n
               DelG(i) = G_2(i) - G_1(i)
               DelX(i) = X_2(i) - X_1(i)
           end do
           do i = 1, n
               HDelG(i) = 0.d0
               do j = 1, n
                   hdelg(i) = hdelg(i) + HI_1(i,j) * DelG(j)
               end do
           end do
           fac = 0.d0
           fae = 0.d0
           sumdg = 0.d0
           sumdx = 0.d0
           do i = 1, n
               fac = fac + delg(i) * delx(i)
               fae = fae + delg(i) * hdelg(i)
               sumdg = sumdg + delg(i)**2
               sumdx = sumdx + delx(i)**2
           end do
           fac = 1.d0 / fac
           fad = 1.d0 / fae
           do i = 1, n
               w(i) = fac * DelX(i) - fad * HDelG(i)
           end do
           DO I = 1, N
               do j = 1, n
                   HI_2(i,j) = HI_1(i,j) + fac * delx(i) * delx(j) -
C change fmc
C&                     fad * HDelG(I) * HDelG(j) + fae * w(i) * w(j)
     1                 fad * HDelG(I) * HDelG(j) + fae * w(i) * w(j)
               end do
           end do
           do i = 1, n
               ChgeX(i) = 0.d0
               do j = 1, n
                   ChgeX(i) = ChgeX(i) - HI_2(i,j) * G_2(j)
               end do
           end do
       END IF

       fac = 0.d0
       do i = 1, n
           fac = fac + ChgeX(i)**2
       end do
       stpl = SQRT(fac)
       IF (stpl .gt. stpmax) THEN
           do i = 1, n
               ChgeX(i) = ChgeX(i) / stpl * stpmax
           end do
       END IF
       lgstst = 0.d0
       do i = 1, n
            IF (ABS(ChgeX(i)) .gt. lgstst) lgstst = ABS(ChgeX(i))
       end do
       IF (lgstst .gt. STPMX) THEN
           do i = 1, n
               ChgeX(i) = ChgeX(i) / lgstst * STPMX
           end do
       END IF
       do i = 1, n
           X_3(i) = X_2(i) + ChgeX(i)
       end do

       END
      SUBROUTINE WriteGeomFile(Natom,Nx,AtNum,X)
      implicit none

      INTEGER i, j, Natom, Nx, AtNum(Natom)
      DOUBLE PRECISION X(Nx)

      OPEN (UNIT=8,FILE="geom")
      DO I = 1, Natom
          write (8,'(I4,3F14.8)') AtNum(I), (X(J),J=3*(I-1)+1,3*(I-1)+3)
      END DO
      write (8,*)
      CLOSE(8)
      return

      END

      SUBROUTINE WriteProgFile(natom,nx,AtNum,Nstep,X_2,X_3,HI,Ea,Eb,Ga,Gb,G)
      implicit none

      INTEGER natom, nx, AtNum(Natom), Nstep
      DOUBLE PRECISION X_2(Nx), X_3(Nx), HI(Nx,Nx), Ea, Eb, Ga(Nx), Gb(Nx), G(Nx)

      INTEGER I, J

      OPEN(UNIT=8,FILE="ProgFile")
      WRITE(8,*) "Progress File for MECP Optimization"
      WRITE(8,*) "Number of Atoms:"
      WRITE(8,*) Natom
      WRITE(8,*) "Number of Steps already Run"
      WRITE(8,*) Nstep
      WRITE(8,*) "Is this a full ProgFile ?"
      WRITE(8,*) 1

      WRITE(8,*) "Next Geometry to Compute:"
      DO I = 1, natom
          WRITE(8,'(I3,3F20.12)') AtNum(I), (X_3(J),J=3*(I-1)+1,3*(I-1)+3)
      END DO

      WRITE(8,*) "Previous Geometry:"
      DO I = 1, Natom
          write (8,'(3F20.12)') (X_2(J),J=3*(I-1)+1,3*(I-1)+3)
      END DO
      WRITE(8,*) "Energies of First, Second State at that Geometry:"
      WRITE(8,'(F20.12)') Ea
      WRITE(8,'(F20.12)') Eb
      WRITE(8,*) "Gradient of First State at that Geometry:"
      DO I = 1, nx
          WRITE(8,'(F20.12)') Ga(I)
      END DO
      WRITE(8,*) "Gradient of Second State at that Geometry:"
      DO I = 1, nx
          WRITE(8,'(F20.12)') Gb(I)
      END DO
      WRITE(8,*) "Effective Gradient at that Geometry:"
      DO I = 1, nx
          WRITE(8,'(F20.12)') G(I)
      END DO
      WRITE(8,*) "Approximate Inverse Hessian at that Geometry:"
      DO I = 1, nx
          DO J = 1, nx
              WRITE(8,'(F20.12)') HI(I,J)
          END DO
      END DO

      CLOSE(8)
      return

      END
"""


if __name__ == '__main__':
    main()
