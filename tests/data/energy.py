def parse_energy(f, line, fields, default=None):
    if len(fields) > 4 and fields[0] == 'SCF' and fields[1] == 'Done:':
        return float(fields[4])
    return default