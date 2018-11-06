#!/usr/bin/env python

from __future__ import print_function
import io
import os
import shutil
import sys
import re
import math
from subprocess import check_output
import pytest
import numpy as np
from easymecp.easymecp import MECPCalculation, temporary_directory
here = os.path.abspath(os.path.dirname(__file__))
data = os.path.join(here, 'data')


AVAILABLE_MEMORY_GB = float(check_output(['free']).splitlines()[1].split()[1]) / 1048576


def required_steps(path):
    step = 0
    with open(path) as f:
        for line in f:
            if 'Geometry at Step' in line:
                step = int(line.split()[-1])
    return step


def too_much_memory(path):
    with open(path) as f:
        for line in f:
            mem = re.search(r'^%mem=([\d\.]+)\S+', line)
            if mem and float(mem.group(1)) > AVAILABLE_MEMORY_GB:
                return True


@pytest.mark.parametrize("directory", sorted(next(os.walk(data))[1]))
def test_easymecp(directory):
    cwd = os.getcwd()
    original_data = os.path.join(here, 'data', directory)

    if 'singlefile' in directory or 'freq' in directory:
        pytest.skip('This entry is meant to be tested by other function')
        return
    if too_much_memory(os.path.join(original_data, 'Input_Header_A')):
        pytest.skip('This entry requires too much memory')
        return
    steps = required_steps(os.path.join(original_data, 'ReportFile'))
    if 'QUICKTEST' in os.environ and steps > 20:
        pytest.skip('Too many steps are expected. Skipping...')
        return

    with temporary_directory() as tmp:
        new_data = os.path.join(tmp, directory)
        shutil.copytree(original_data, new_data)
        os.chdir(new_data)
        kwargs = {'geom': 'geom_init', 'footer': 'Input_Tail'}
        if directory == 'C6H5+_B1-3A2':
            kwargs['energy_parser'] = 'td'
        calc = MECPCalculation(**kwargs)
        result = calc.run()
        if result != calc.OK:
            try:
                os.makedirs(os.path.join(cwd, 'failed'))
            except:
                pass
            shutil.copytree(new_data, os.path.join(cwd, 'failed', directory))
        assert result == calc.OK
        if calc.converged_at != steps:
            print('! Warning: calculation ended OK but took a different number of steps')


def test_external_energy():
    cwd = os.getcwd()
    directory = 'CH2'
    with temporary_directory() as tmp:
        original_data = os.path.join(here, 'data', directory)
        energy = os.path.join(here, 'data', 'energy.py')
        steps = required_steps(os.path.join(original_data, 'ReportFile'))
        new_data = os.path.join(tmp, directory)
        shutil.copytree(original_data, new_data)
        os.chdir(new_data)
        calc = MECPCalculation(geom='geom_init', footer='Input_Tail', energy_parser=energy)
        result = calc.run()
        if result != calc.OK:
            try:
                os.makedirs(os.path.join(cwd, 'failed'))
            except:
                pass
            shutil.copytree(new_data, os.path.join(cwd, 'failed', directory))
        assert result == calc.OK
        if calc.converged_at != steps:
            print('! Warning: calculation ended OK but took a different number of steps')


def test_singlefile():
    directory = 'C6H5+_singlefile'
    cwd = os.getcwd()
    original_data = os.path.join(here, 'data', directory)
    steps = required_steps(os.path.join(original_data, 'ReportFile'))
    with temporary_directory() as tmp:
        new_data = os.path.join(tmp, directory)
        shutil.copytree(original_data, new_data)
        os.chdir(new_data)
        calc = MECPCalculation.from_gaussian_input_file('input.gjf')
        result = calc.run()
        if result != calc.OK:
            try:
                os.makedirs(os.path.join(cwd, 'failed'))
            except:
                pass
            shutil.copytree(new_data, os.path.join(cwd, 'failed', directory))
        assert result == calc.OK
        if calc.converged_at != steps:
            print('! Warning: calculation ended OK but took a different number of steps')


@pytest.mark.parametrize("directory, freq_a, freq_b, energy_a, energy_b, energy_avg", [
    ('C6H5+_freq', -552.2857, 317.7488, -231.188991, -231.186626, -231.1878085),
])
def test_freq(directory, freq_a, freq_b, energy_a, energy_b, energy_avg):
    cwd = os.getcwd()
    original_data = os.path.join(here, 'data', directory)
    steps = required_steps(os.path.join(original_data, 'ReportFile'))
    with temporary_directory() as tmp:
        new_data = os.path.join(tmp, directory)
        shutil.copytree(original_data, new_data)
        os.chdir(new_data)
        calc = MECPCalculation.from_gaussian_input_file('input.gjf')
        result = calc.run()
        if result != calc.OK:
            try:
                os.makedirs(os.path.join(cwd, 'failed'))
            except:
                pass
            shutil.copytree(new_data, os.path.join(cwd, 'failed', directory))
        assert result == calc.OK
        if calc.converged_at != steps:
            print('! Warning: calculation ended OK but took a different number of steps')

        with open('ReportFile') as f:
            for line in f:
                if line.startswith('Minimum frequency for state A'):
                    value = float(line.split()[-2])
                    assert abs(freq_a - value) < 1e-3
                elif line.startswith('Minimum frequency for state B'):
                    value = float(line.split()[-2])
                    assert abs(freq_b - value) < 1e-3
                elif line.startswith('Sum of electronic and thermal Free Energies for state B'):
                    value = float(line.split()[-2])
                    assert abs(energy_a - value) < 1e-2
                elif line.startswith('Sum of electronic and thermal Free Energies for state A'):
                    value = float(line.split()[-2])
                    assert abs(energy_b - value) < 1e-2
                elif line.startswith('Avg sum of electronic and thermal Free Energies for both states'):
                    value = float(line.split()[-2])
                    assert abs(energy_avg - value) < 1e-2


def test_geometry():
    directory = 'C6H5+'
    # directory = 'C6H5+_singlefile'

    # Values provided in Theor Chem Acc (1998) 99:95-99, Table 2, MECP B3LYP row
    distances = [1.415, 1.392, 1.437]
    angles = [128.9, 114.9, 119.4, 122.5, 122.3, 120.6]

    cwd = os.getcwd()
    original_data = os.path.join(here, 'data', directory)
    steps = required_steps(os.path.join(original_data, 'ReportFile'))
    with temporary_directory(remove=False) as tmp:
        new_data = os.path.join(tmp, directory)
        print(new_data, file=sys.stderr)
        shutil.copytree(original_data, new_data)
        os.chdir(new_data)
        calc = MECPCalculation(geom='geom_init')
        # calc = MECPCalculation.from_gaussian_input_file('input.gjf')
        result = calc.run()
        if result != calc.OK:
            try:
                os.makedirs(os.path.join(cwd, 'failed'))
            except:
                pass
            shutil.copytree(new_data, os.path.join(cwd, 'failed', directory))
        assert result == calc.OK
        if calc.converged_at != steps:
            print('! Warning: calculation ended OK but took a different number of steps')

        if os.path.isfile('geom'):
            atoms = parse_xyz('geom')
        else:
            atoms = parse_xyz(calc.geom)
        # List contains: C1, C2, C3, C4, C5, C6, H1, H2, H3, H4, H5; in that order
        assert abs(distances[0] - np_distance(atoms[0][1], atoms[5][1])) < 0.15
        assert abs(distances[1] - np_distance(atoms[4][1], atoms[5][1])) < 0.15
        assert abs(distances[2] - np_distance(atoms[3][1], atoms[4][1])) < 0.15
        assert abs(angles[0] - angle(atoms[5][1], atoms[0][1], atoms[1][1])) < 1.5
        assert abs(angles[1] - angle(atoms[0][1], atoms[1][1], atoms[2][1])) < 1.5
        assert abs(angles[2] - angle(atoms[1][1], atoms[2][1], atoms[3][1])) < 1.5
        assert abs(angles[3] - angle(atoms[2][1], atoms[3][1], atoms[4][1])) < 1.5
        assert abs(angles[4] - angle(atoms[0][1], atoms[1][1], atoms[6][1])) < 1.5
        assert abs(angles[5] - angle(atoms[1][1], atoms[2][1], atoms[7][1])) < 1.5
        assert rmsd('geom_init', calc.geom) < 0.1
        with open(os.path.join(cwd, 'geometry-results.txt'), 'w') as f:
            print('Distance', atoms[0][0], atoms[5][0], '=', np_distance(atoms[0][1], atoms[5][1]), file=f)
            print('Distance', atoms[4][0], atoms[5][0], '=', np_distance(atoms[4][1], atoms[5][1]), file=f)
            print('Distance', atoms[3][0], atoms[4][0], '=', np_distance(atoms[3][1], atoms[4][1]), file=f)
            print('Angle', atoms[5][0], atoms[0][0], atoms[1][0], '=', angle(atoms[5][1], atoms[0][1], atoms[1][1]), file=f)
            print('Angle', atoms[0][0], atoms[1][0], atoms[2][0], '=', angle(atoms[0][1], atoms[1][1], atoms[2][1]), file=f)
            print('Angle', atoms[1][0], atoms[2][0], atoms[3][0], '=', angle(atoms[1][1], atoms[2][1], atoms[3][1]), file=f)
            print('Angle', atoms[2][0], atoms[3][0], atoms[4][0], '=', angle(atoms[2][1], atoms[3][1], atoms[4][1]), file=f)
            print('Angle', atoms[0][0], atoms[1][0], atoms[6][0], '=', angle(atoms[0][1], atoms[1][1], atoms[6][1]), file=f)
            print('Angle', atoms[1][0], atoms[2][0], atoms[7][0], '=', angle(atoms[1][1], atoms[2][1], atoms[7][1]), file=f)
            print('RMSD =', rmsd('geom_init', calc.geom), file=f)


def parse_xyz(path):
    atoms = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if len(fields) == 4:
                atoms.append((fields[0], (float(fields[1]), float(fields[2]), float(fields[3]))))
    return atoms


def distance(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def np_distance(a, b):
    return np.linalg.norm(np.array(a)-np.array(b))


def angle(a, b, c):
    # Create vectors from points
    ba = [aa-bb for (aa,bb) in zip(a,b)]
    bc = [cc-bb for (cc,bb) in zip(c,b)]

    # Normalize vector
    nba = math.sqrt(sum((x**2.0 for x in ba)))
    ba = [x/nba for x in ba]

    nbc = math.sqrt(sum((x**2.0 for x in bc)))
    bc = [x/nbc for x in bc]

    # Calculate scalar from normalized vectors
    scal = sum((aa*bb for (aa,bb) in zip(ba,bc)))

    # calculate the angle in radian
    angle = math.acos(scal)
    # return in degrees
    return angle * 180 / math.pi


def rmsd(file_a, file_b):
    mol_a = parse_xyz(file_a)
    mol_b = parse_xyz(file_b)

    xyz_a = np.array([atom[1] for atom in mol_a])
    xyz_b = np.array([atom[1] for atom in mol_b])

    return np.sqrt(np.mean(np.square(xyz_a-xyz_b)))