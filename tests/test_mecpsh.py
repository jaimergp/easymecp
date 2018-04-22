#!/usr/bin/env python

import os
import shutil
import pytest
from mecpsh.mecpsh import MECPCalculation, temporary_directory
here = os.path.abspath(os.path.dirname(__file__))

@pytest.mark.parametrize("directory",
sorted(os.listdir(os.path.join(here, 'data')))
)
def test_mecpsh(directory):
    cwd = os.getcwd()
    with temporary_directory() as tmp:
        shutil.copytree(os.path.join(here, 'data', directory), os.path.join(tmp, directory))
        os.chdir(directory)
        calc = MECPCalculation(geom='geom_init')
        result = calc.run()
        assert result == calc.OK