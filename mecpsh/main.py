#!/usr/bin/env python

"""
MECPsh
------

Simplified Python wrapper around the original MECP Fortran code
from J. Harvey (2003).

It uses the same source code, but saves the hassle of manual
edition of files and shell loops. Now, you only* need to provide:

- header contents (formerly Input_Header_A and Input_Header_B)
- starting geometry (geom file)
- footer contents (footer file)

* number of atoms is inferred from geometry file automatically

These options must be specified in an inputfile formatted like this:

    a_header: filename_a
    b_header: filename_b
    ... etc

... or directly with commandline arguments. Available configuration
keys can be always printed with `-h` flag.
"""

from __future__ import print_function, absolute_import
from datetime import datetime
from subprocess import call
from tempfile import TemporaryDirectory
import argparse
import os
import shutil
import sys


class MECPCalculation(object):

    OK = 'OK'
    ERROR = 'ERROR'
    MAX_ITERATIONS_REACHED = 'MAX_ITERATIONS_REACHED'

    def __init__(self, max_steps=20, a_header='Input_Header_A', b_header='Input_Header_B',
                 geom='geom', footer='footer', gaussian_exe='g09', TDE='5.d-5',
                 TDXMax='4.d-3', TDXRMS='2.5d-3', TGMax='7.d-4', TGRMS='5.d-4',
                 FC='gfortran', FFLAGS='-O -ffixed-line-length-none', energy_parser='dft', **kwargs):
        self.max_steps = max_steps
        self.a_header = a_header
        self.b_header = b_header
        self.geom = geom
        self.footer = footer
        self.gaussian_exe = gaussian_exe
        self.TDE = TDE
        self.TDXMax = TDXMax
        self.TDXRMS = TDXRMS
        self.TGMax = TGMax
        self.TGRMS = TGRMS
        self.FC = FC
        self.FFLAGS = FFLAGS
        self.natom = 0
        if energy_parser in ('dft', 'mp2'):
            self._parse_energy = getattr(self, '_parse_energy_' + energy_parser)
        else:
            raise ValueError('energy_parser must be `dft` or `mp2`')
        with open(geom) as f:
            for line in f:
                if len(line.split()) == 4:
                    self.natom +=1
    @classmethod
    def from_conf(cls, path):
        d = _get_defaults()
        with open(path) as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#') or ':' not in line:
                    continue
                fields = line.split(':')
                if len(fields) < 2:
                    print('Malformed line #{}: {}'.format(i, line))
                    sys.exit()
                key, value = fields[0], ':'.join(fields[1:])
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

    def run(self):
        # Recompile Fortran code
        print('Preparing workspace...')
        self.prepare_workspace()
        print('Compiling MECP...')
        mecp_exe = self.recompile_fortran()
        geom = self.geom

        for i in range(self.max_steps):
            # First, compute and parse Gaussian Jobs
            print('Running step #{}...'.format(i))
            print('Launching Gaussian job for file A...')
            input_a = self.prepare_gaussian(self.a_header, geom, self.footer, label='A', step=i)
            energy_a, gradients_a = self.run_gaussian(input_a)
            print('Launching Gaussian job for file B...')
            input_b = self.prepare_gaussian(self.b_header, geom, self.footer, label='B', step=i)
            energy_b, gradients_b = self.run_gaussian(input_b)

            # Second, run MECP
            if not self.prepare_ab_initio(energy_a, energy_b, gradients_a, gradients_b):
                return self.ERROR
            print('Launching MECP...')
            try:
                call([mecp_exe])
            except Exception as e:
                print('  ! Error during MECP execution:', type(e), e)
                self.report('ERROR')
            if os.path.isfile('AddtoReportFile'):
                with open('AddtoReportFile') as f:
                    self.report(f.read())
                os.remove('AddtoReportFile')

            # Third, check results in this iteration
            with open('ReportFile') as f:
                for line in f:
                    if 'CONVERGED' in line:
                        print('MECP optimization has converged at Step', i)
                        return self.OK
                    if 'ERROR' in line:
                        print('An error ocurred!')
                        return self.ERROR
            # Keep going! Next iteration will use geom autogenerated by MECP.x
            geom = 'geom'

        return self.MAX_ITERATIONS_REACHED

    def recompile_fortran(self):
        with TemporaryDirectory() as tmpdir:
            tmpsource = os.path.join(tmpdir, 'src')
            shutil.copytree(_SOURCE, tmpsource)
            # Patch number of atoms
            with open(os.path.join(tmpsource, 'MainProgram.f'), 'r+') as f:
                contents = f.read().replace('parameter (Natom=77)',
                                            'parameter (Natom={})'.format(self.natom))
                f.seek(0)
                f.write(contents)
                f.truncate()
            # Patch convergence criteria
            with open(os.path.join(tmpsource, 'TestConvergence.f'), 'r+') as f:
                search = 'PARAMETER (TDE=5.d-5,TDXMax=4.d-3,TDXRMS=2.5d-3,TGMax=7.d-4,TGRMS=5.d-4))'
                replacement = 'PARAMETER (TDE={},TDXMax={},TDXRMS={},TGMax={},TGRMS={}))'.format(
                    self.TDE, self.TDXMax, self.TDXRMS, self.TGMax, self.TGRMS)
                contents = f.read().replace(search, replacement)
                f.seek(0)
                f.write(contents)
                f.truncate()
            # Patch compiler options
            with open(os.path.join(tmpsource, 'makefile'), 'r+') as f:
                contents = f.read().replace('FC = gfortran', 'FC = {}'.format(self.FC))
                contents = contents.replace('FFLAGS= -O -ffixed-line-length-none', 'FFLAGS= {}'.format(self.FFLAGS))
                f.seek(0)
                f.write(contents)
                f.truncate()
            with open('__fortran_compilation', 'w') as out:
                call(['make'], cwd=tmpsource, stdout=out)
            shutil.copyfile(os.path.join(tmpsource, 'MECP.x'), 'MECP.x')
            os.chmod('MECP.x', os.stat('MECP.x').st_mode | 0o111)  # make executable
        return './MECP.x'

    def element_symbol_to_number(self, fh):
        elements = ELEMENTS
        lines = []
        for line in fh:
            fields = line.split()
            if len(fields) > 3:
                if not fields[0].isdigit():
                    number = elements.get(fields[0].title(), -1)
                    line = line.replace(fields[0], str(number))
            lines.append(line)
        return '\n'.join(lines)

    def prepare_workspace(self):
        i = 1
        while os.path.exists('JOBS'):
            try:
                os.rename('JOBS', 'JOBS{}'.format(i))
            except:
                i += 1
        os.makedirs('JOBS')
        with open(self.geom) as f:
            geometry = self.element_symbol_to_number(f)
        with open('ProgFile', 'w') as f:
            f.write(_PROGFILE.format(natom=self.natom, geometry=geometry))
        with open('ReportFile', 'w') as f:
            print(datetime.now(), file=f)

    def prepare_gaussian(self, header, geom, footer, label='A', step=0):
        name = 'Job{}_{}.gjf'.format(step, label)
        with open(name, 'w') as f:
            with open(header) as a:
                f.write(a.read())
            with open(geom) as b:
                contents = self.element_symbol_to_number(b)
                f.write(contents)
            with open(footer) as c:
                f.write(c.read())
        return name

    def run_gaussian(self, inputfile):
        try:
            call([self.gaussian_exe, inputfile])
        except:
            print('  ! Could not run Gaussian job', inputfile)
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
        if fields[0] == 'Center' and fields[1] == 'Atomic' and fields[2] == 'Forces':
            gradients = []
            next(f), next(f)  # skip two lines
            for i in range(self.natom):
                line = next(f)
                gradients.append(line.split()[1:5])
            return gradients
        return default

    def _parse_energy_dft(self, f, line, fields, default=None):
        if fields[0] == 'SCF' and fields[1] == 'Done:':
            return float(fields[4])
        return default

    def _parse_energy_mp2(self, f, line, fields, default=None):
        if fields[0] == 'E2' and fields[3] == 'EUMP2':
            return float(fields[5])
        return default

    def prepare_ab_initio(self, energy_a, energy_b, gradients_a, gradients_b):
        if not all(x is not None for x in (energy_a, energy_b, gradients_a, gradients_b)):
            print('  ! Some energies or gradients could not be obtained!')
            return False
        with open('ab_initio', 'w') as f:
            print('Energy of the First State', file=f)
            print(energy_a, file=f)
            print('Gradient of the First State' , file=f)
            print('\n'.join('{}  {}   {}   {}'.format(*[l for l in gradients_a])), file=f)
            print('Energy of the Second State' , file=f)
            print(energy_b, file=f)
            print('Gradient of the Second State' , file=f)
            print('\n'.join('{}  {}   {}   {}'.format(*[l for l in gradients_b])), file=f)
        return True

    def report(self, msg):
        with open('ReportFile', 'a') as r:
            print(msg, file=r)


def _parse_cli():
    p = argparse.ArgumentParser(prog='mecpsh', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-f', '--inputfile', metavar='FILE', help='Configuration file.')
    defaults = _get_defaults()
    for k, v in sorted(defaults.items()):
        p.add_argument('--'+k, default=v, metavar='VALUE', help='{} (default={!r})'.format(USAGE[k], v))
    args = p.parse_args()
    return args


def _valid_fortran_double(value, key=None):
    try:
        float(value.replace('d', 'e'))
    except ValueError:
        raise ValueError('Value {}{} must be a valid Fortran double ([base].d[exp], like 5.d-4)'.format(
                         value, 'for {}'.format(key) if key is not None else ''))
    else:
        return value


def main():
    args = _parse_cli()
    if args.inputfile:
        print('Configuration file', args.inputfile, 'has been specified. Ignoring any other arguments!')
        calc = MECPCalculation.from_conf(args.inputfile)
    else:
        calc = MECPCalculation(**args.__dict__)
    result = calc.run()
    if result == MECPCalculation.OK:
        print('Success!')
    elif result == MECPCalculation.ERROR:
        print('Something failed...')
    elif result == MECPCalculation.MAX_ITERATIONS_REACHED:
        print('Calculation did not converge after', calc.max_steps, 'steps...')


def _get_defaults():
    _defaults = MECPCalculation.__init__.__defaults__
    _ndef = len(_defaults)
    _narg = MECPCalculation.__init__.__code__.co_argcount
    _args = MECPCalculation.__init__.__code__.co_varnames[(_narg-_ndef):_narg]
    return {k:v for (k, v) in zip(_args, _defaults)}


DEFAULTS = _get_defaults()
_SOURCE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src')
_PROGFILE = """
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
USAGE = {
    'max_steps': 'Number of iterations to perform until stop or convergence',
    'a_header': 'File containing the top part of system configuration with multiplicity A',
    'b_header': 'File containing the top part of system configuration with multiplicity B',
    'geom': 'File containing the starting system geometry (element symbols will be converted in atomic numbers automatically)',
    'footer': 'File containing the bottom part of system configuration',
    'energy_parser': 'Which energy should be parsed: dft or mp2.',
    'gaussian_exe': 'Path to gaussian executable. Compatible versions: g09, g16',
    'TDE': 'Convergence criterion. Must be a valid Fortran double!',
    'TDXMax': 'Convergence criterion. Must be a valid Fortran double!',
    'TDXRMS': 'Convergence criterion. Must be a valid Fortran double!',
    'TGMax': 'Convergence criterion. Must be a valid Fortran double!',
    'TGRMS': 'Convergence criterion. Must be a valid Fortran double!',
    'FC': 'Fortran compiler',
    'FFLAGS': 'Fortran compiler flags'
}

if __name__ == '__main__':
    main()