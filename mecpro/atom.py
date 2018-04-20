#!/usr/bin/python

#Copyright (c) 2015 Brigham Young University

#See the file license.txt for copying permission.

import itertools
import math

import rotate

from copy import deepcopy

elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf',
'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Rg', 'Hs']

atomic_numbers = {}

i = 1
for element in elements:
    atomic_numbers[element] = i
    i += 1

'''Representation of an atom. Has an element and coordinates'''
class Atom (object):

    def __init__(self, element, x, y, z):
        self.x, self.y, self.z = float(x), float(y), float(z)
        if isinstance(element, int) or element.isdigit():
            self.element = elements[int(element) - 1]
        else:
            assert(element in elements)
            self.element =  element

    @property
    def coords(self):
        return self.x, self.y, self.z

    @coords.setter
    def coords(self, value):
        self.x, self.y, self.z = tuple(value)

    def __str__(self):
        return "%(element)s\t%(x) .8f\t%(y) .8f\t%(z) .8f" % self.__dict__

    def __repr__(self):
        return "<Atom (%(element)s: %(x)f, %(y)f, %(z)f)>" % self.__dict__

'''Molecules are lists of atoms with some additional operations'''
class Molecule (object):

    def __init__(self, atoms = None):
        if atoms is None:
            self.atoms = []
        elif isinstance(atoms, Molecule):
            self.atoms = deepcopy(atoms.atoms)
        else:
            self.atoms = deepcopy(atoms)

    def append(self, atom):
        self.atoms.append(atom)

    @property
    def coords(self):
        return map(lambda x: x.coords, self.atoms)

    @coords.setter
    def coords(self, value):
        for atom, value in itertools.izip(self.atoms, value):
            atom.coords = value

    '''Returns a list of all the coordinates'''
    @property
    def components(self):
        return list(itertools.chain(*self.coords))

    @components.setter
    def components(self, value):
        coords = []
        current = []
        count = 0
        for comp in value:
            current.append(comp)
            count += 1
            if count == 3:
                coords.append(tuple(current))
                count = 0
                current = []
        self.coords = coords

    def __str__(self):
        return '\n'.join(map(str, self.atoms))

    def __repr__(self):
        return '<Molecule (' + ', '.join(map(repr, self.atoms)) + ')>'

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def __getitem__(self, index):
        return self.atoms[index]

    def __setitem__(self, index, value):
        if not isinstance(value, Atom):
            raise TypeError("Not an atom: " + value)
        self.atoms[index] = value

    '''Normalize the geometry of this molecule:
Translate and rotate it such that:
* The first atom is at the origin
* The second atom is on the y-axis
* The third atom is on the x-y plane'''
    def normalize(self):
        coords = self.coords
        trans = map(lambda x: -x, coords[0])
        # translate everything so that atom 1 is at (0,0,0)
        for i in xrange(len(coords)):
            coords[i] = rotate.add_vec(coords[i], trans)

        if len(coords) < 2:
            self.coords = coords
            return

        # Rotate everything so that atom 2 is on the y axis
        axis = rotate.get_rot_axis((0,0,0), coords[1], (0,1,0))
        angle = rotate.get_rot_angle(coords[1], (0,1,0))
        coords = rotate.rotate_many(coords, angle, (0,0,0), axis)

        if len(coords) < 3:
            self.coords = coords
            return

        # Rotate everything so that atom 3 is on the x-y plane
        axis = (0,1,0) # rotate around the y axis, where the existing 2 normalized points are
        angle = math.atan2(-coords[2][2], coords[2][0]) # (x, -z)
        self.coords = rotate.rotate_many(coords, angle, (0,0,0), axis)

# small unit test
# _mol_test = Molecule([Atom('H', 1.0, 0.0, 2.0), Atom('H', -1.0, 0.0, 2.0),
                      # Atom('O', 0.0, 1.0, 2.1), Atom('C', 1.3, 1.5, 2.3)])
# _mol_test.normalize()
# print _mol_test

# assert _mol_test.coords == [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
# _mol_test.coords         = [(1.0, 0.0, 0.0), (-1.0, 2.0, 0.0), (0.0, 1.0, 0.0)]
# assert _mol_test.coords == [(1.0, 0.0, 0.0), (-1.0, 2.0, 0.0), (0.0, 1.0, 0.0)]
# assert _mol_test.atoms[1].y == 2.0

# assert _mol_test.components == [1.0, 0.0, 0.0, -1.0, 2.0, 0.0, 0.0, 1.0, 0.0]
# _mol_test.components         = [1.0, 0.5, 0.0, -1.0, 2.0, 0.0, 0.0, 1.0, 0.0]
# assert _mol_test.components == [1.0, 0.5, 0.0, -1.0, 2.0, 0.0, 0.0, 1.0, 0.0]
# assert _mol_test.coords == [(1.0, 0.5, 0.0), (-1.0, 2.0, 0.0), (0.0, 1.0, 0.0)]

#del _mol_test
