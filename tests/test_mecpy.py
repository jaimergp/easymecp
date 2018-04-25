#!/usr/bin/env python

import io
import os
import shutil
import sys
import pytest
from mecpy.mecpy import MECPCalculation, temporary_directory
here = os.path.abspath(os.path.dirname(__file__))


def required_steps(path):
    with open(path) as f:
        for line in f:
            if 'Geometry at Step' in line:
                step = int(line.split()[-1])
    return step


@pytest.mark.parametrize("directory", sorted(next(os.walk(os.path.join(here, 'data')))[1]))
def test_mecpsh(directory):
    cwd = os.getcwd()
    original_data = os.path.join(here, 'data', directory)
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