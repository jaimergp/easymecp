#!/usr/bin/python

#Copyright (c) 2015 Brigham Young University

#See the file license.txt for copying permission.

# Needed to interface with Gaussian
import subprocess

import re
import atom

from datetime import datetime

_gaussian = "/fslapps/chem/bin/rung09"

_gprocess = None

_blank = lambda line: (line == '' or line.isspace())

def linestrip(lst):
    while len(lst) > 0 and _blank(lst[0]):
        del lst[0]

    while len(lst) > 0 and _blank(lst[-1]):
        del lst[-1]

    return [s.rstrip() for s in lst]

'''If an instance of gaussian is currently running, kill it'''
def kill():
    if _gprocess is not None:
        _gprocess.kill()

'''If an instance of gaussian is currently running, terminate it'''
def terminate():
    if _gprocess is not None:
        _gprocess.terminate()

'''Generic generator for making header'''
def _header_template(title, charge, spin, type = '', basis = 'm06/gen', command = {}, link0 = {}):
    options = []
    for k, v in link0.iteritems():
        yield '%' + k + '=' + str(v) + '\n'

    options = []
    for k, v in command.iteritems():
        if v is None:
            options.append(k)
        elif not isinstance(v, (list, tuple)):
            options.append(k + '=' + v)
        else:
            v = map(str, v)
            options.append(k + '(' + ','.join(v) + ')')

    yield '#%s %s %s\n' % (type, basis, ' '.join(options))

    yield '\n'
    yield title + '\n'
    yield '\n'

    yield '%d %d\n' % (charge, spin)

'''Creates the input files for gaussian to run.'''
def write_input(name, molecule, header, footer = []):
    with open(name + ".com", 'w') as f:
        f.write('\n'.join(linestrip([l.rstrip() for l in header])))
        f.write('\n')
        f.write(str(molecule))
        f.write('\n\n')
        tail = linestrip([l.rstrip() for l in footer])
        if len(tail) > 0:
            f.write('\n'.join(tail))
            f.write('\n\n')

_geom_marker = re.compile(r"\d\s+\d")
_float_pat = r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
_geom_pat = re.compile(r'\s*(?:' + '|'.join(atom.elements) + r'|\d+)\s+' +
    r'\s+'.join([_float_pat] * 3) + r'\s*')

def read_input(name):
    meta = {'header': [], 'footer': []} # before and after the geometries
    molecule = atom.Molecule()
    part = 'header'
    with open(name + '.com') as f:
        for line in f:

            if part == "molecule":
                if _geom_pat.match(line):
                    element, x, y, z = line.split()
                    coords = map(float, (x, y, z))
                    molecule.append(atom.Atom(element, *coords))
                elif _blank(line) and len(molecule) == 0:
                    continue
                else:
                    part = 'footer'
                    meta[part].append(line)
            else:
                meta[part].append(line)
                if _geom_marker.match(line):
                    part = "molecule"

    return molecule, meta['header'], meta['footer']

"""Create a Gaussian input file and run it.
Returns the status code and the running time. (Running time to be used for profiling)"""
def run_g09(name, title, charge, spin, molecule,
            type = '', basis = 'm06/gen', command = {},
            footer = [], link0 = {}):
    header = _header_template(title, charge, spin, type, basis, command, link0)
    write_input(name, molecule, header, footer)
    start_time = datetime.now()
    global _gprocess
    _gprocess = subprocess.Popen([_gaussian, name + '.com'])
    status = _gprocess.wait()
    _gprocess = None
    return status, datetime.now() - start_time
