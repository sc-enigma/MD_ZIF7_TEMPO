import numpy as np

def read_mol2_file(filename):
    mol2file = open(filename, 'r')
    mol2sections = {}

    current_section = ''
    current_lines = []
    for line in mol2file:
        if '@' in line:
            if len(current_lines) > 0:
                mol2sections[current_section] = current_lines.copy()
            current_lines.clear()
            current_section = line.replace('@<TRIPOS>', '').replace('\n', '')
        else:
            current_lines.append(line)
    mol2file.close()
    return mol2sections