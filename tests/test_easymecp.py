#!/usr/bin/env python

import io
import os
import shutil
import sys
import re
from subprocess import check_output
import pytest
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