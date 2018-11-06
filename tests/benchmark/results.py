import numpy as np

def parse_xyz(path):
    atoms = None
    with open(path) as f:
        for line in f:
            if line.startswith('Geometry at'):
                atoms = []
                line = next(f)
                while line.strip():
                    fields = line.split()
                    if len(fields) == 4:
                        atoms.append((fields[0], [float(i) for i in fields[1:]]))
                    line = next(f)
    return atoms

def rmsd(file_a, file_b):
    mol_a = parse_xyz(file_a)
    mol_b = parse_xyz(file_b)

    xyz_a = np.array([atom[1] for atom in mol_a])
    xyz_b = np.array([atom[1] for atom in mol_b])

    return np.sqrt(np.mean(np.square(xyz_a-xyz_b)))


def energy(path):
    energy = None
    with open(path) as f:
        for line in f:
            if 'Energy of First State' in line:
                energy = float(line.split()[-1])
    return energy


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print('Usage: results.py ReportFileA ReportFileB')
        sys.exit()
    print('ENERGY A ({}) ='.format(sys.argv[1]), energy(sys.argv[1]))
    print('ENERGY B ({}) ='.format(sys.argv[2]), energy(sys.argv[2]))
    print('ENERGY A-B =', abs(energy(sys.argv[1]) - energy(sys.argv[2])))
    print('RMSD =', rmsd(sys.argv[1], sys.argv[2]), 'A')

