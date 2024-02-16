import numpy as np

from atom import Atom
from pbc import calculate_lattice_vectors
from get_topology_params import getAllBonds, getAllAngles, getAllDihedrals

def format_val(val, length_before_comma, length_after_comma=-1, allign_left=False):
    length_required = length_before_comma
    if length_after_comma >= 0:
        val = f'%.{length_after_comma}f' % val
        length_required = length_before_comma + 1 + length_after_comma
    val = str(val)
    if len(val) > length_required:
        val = val[0:length_required]  
    while len(val) < length_required:
        if allign_left:
            val += ' '
        else:
            val = ' ' + val
    return val
        
def write_gro_file(atoms, filename, \
                   a, b, c, \
                   alpha=np.pi/2, beta=np.pi/2, gamma=np.pi/2,
                   bounds_a=[0,1], bounds_b=[0,1], bounds_c=[0,1],
                   scale=1.0):
    file = open(filename, 'w')
    file.write('ITC 2023. ALL RIGHTS RESERVED.\n')
    file.write(f'{format_val(len(atoms), 5)}\n')
    
    # Atoms
    for atom_idx in range(len(atoms)):
        atom = atoms[atom_idx]
        file.write(format_val(atom.resid_idx,          5))
        file.write(format_val(atom.resid,              5, allign_left=True))
        file.write(format_val(atom.name,               5))
        file.write(format_val(atom_idx + 1,            5))
        file.write(format_val(atom.r[0] * 0.1 * scale, 4, 3))
        file.write(format_val(atom.r[1] * 0.1 * scale, 4, 3))
        file.write(format_val(atom.r[2] * 0.1 * scale, 4, 3))
        file.write('\n')
    
    # Box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
    # The last 6 values may be omitted (they will be set to zero)
    # Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
    
    vec_a, vec_b, vec_c = calculate_lattice_vectors(a, b, c, alpha, beta, gamma)
    vec_a, vec_b, vec_c = \
        vec_a * (bounds_a[1] - bounds_a[0]) * 0.1 * scale, \
        vec_b * (bounds_b[1] - bounds_b[0]) * 0.1 * scale, \
        vec_c * (bounds_c[1] - bounds_c[0]) * 0.1 * scale
        # convert to nm
    
    file.write(format_val(vec_a[0], 4, 5))
    file.write(format_val(vec_b[1], 4, 5))
    file.write(format_val(vec_c[2], 4, 5))
    file.write(format_val(0.0,      4, 5))
    file.write(format_val(0.0,      4, 5))
    file.write(format_val(vec_b[0], 4, 5))
    file.write(format_val(0.0,      4, 5))
    file.write(format_val(vec_c[0], 4, 5))
    file.write(format_val(vec_c[1], 4, 5))
    file.write('\n')
    
    file.close()

def write_mol2_file(atoms, filename, a, b, c, alpha, beta, gamma, skip_long_bonds=False):
    bonds = getAllBonds(atoms)
    is_short = np.empty(len(bonds), dtype=bool)
    if skip_long_bonds:
        for bond_idx in range(len(bonds)):
            delta = atoms[bonds[bond_idx][0]].r - atoms[bonds[bond_idx][1]].r
            is_short[bond_idx] = pow(sum(np.power(delta, 2)), 0.5) < 10.0
        bonds = bonds[is_short]
        
    file = open(filename, 'w')
    file.write('@<TRIPOS>MOLECULE\n')
    file.write('MOL\n')
    file.write(f'{format_val(len(atoms), 6)}{format_val(len(bonds), 6)}{format_val(0, 6)}\n')
    file.write('SMALL\n')
    file.write('NO_CHARGES\n')
    file.write('****\n')
    file.write('ITC 2023. ALL RIGHTS RESERVED.\n\n')
    
    file.write('@<TRIPOS>ATOM\n')
    for atom in atoms:
        file.write(format_val(atom.atom_idx+1, 6))
        file.write(' ')
        file.write(format_val(atom.name, 6, allign_left=True))
        file.write(format_val(atom.r[0], 4, 4))
        file.write(format_val(atom.r[1], 4, 4))
        file.write(format_val(atom.r[2], 4, 4))
        file.write('   ')
        file.write(atom.name[0])
        file.write(format_val(atom.resid_idx, 10))
        file.write(' RES')
        file.write(format_val(atom.resid_idx, 4, allign_left=True))
        file.write('0.0000\n')
    file.write('@<TRIPOS>BOND')
    for bond_idx in range(len(bonds)):
        file.write(format_val(bond_idx+1, 6))
        file.write(format_val(bonds[bond_idx][0]+1, 6))
        file.write(format_val(bonds[bond_idx][1]+1, 6))
        file.write('   un\n')
    file.write('@<TRIPOS>CRYSIN\n')
    file.write(format_val(a, 5, 4))
    file.write(format_val(b, 5, 4))
    file.write(format_val(c, 5, 4))
    file.write(format_val(alpha, 5, 4))
    file.write(format_val(beta, 5, 4))
    file.write(format_val(gamma, 5, 4))
    file.write('\n')
    file.close()

def write_atoms(atoms, charge, mass, resid, filename):
    file = open(filename, 'w')
    file.write('[ atoms ]\n')
    file.write(';   nr       type  resnr residue  atom   cgnr     charge       mass\n')
    
    for atom_idx in range(len(atoms)):
        file.write(format_val(atom_idx+1, 6))
        file.write(format_val(atoms[atom_idx].atom_type, 11))
        file.write(format_val(1, 7))
        file.write(format_val(resid, 7))
        file.write(format_val(atoms[atom_idx].name, 6))
        file.write(format_val(1, 7))
        file.write(format_val(charge[atoms[atom_idx].atom_type], 6, 4))
        file.write(format_val(mass[atoms[atom_idx].atom_type], 6, 4))
        file.write('\n')
    file.close()
    
def reverse_param_key(string):
    split = np.flip(string.replace('-', ' ').split())
    string_reversed = split[0]
    for val in split[1:]:
        string_reversed += '-' + val
    return string_reversed

def write_bonds(atoms, bond_params, filename, skip_long_bonds=False):
    bonds = getAllBonds(atoms)
    
    file = open(filename, 'w')
    file.write('[bonds]\n')
    file.write(';   ai     aj funct   r             k\n')
    for bond in bonds:
        if skip_long_bonds:
            delta = np.sqrt(sum(np.power((atoms[bond[0]].r - atoms[bond[1]].r), 2)))
            if delta > 10.0:
                continue
        
        # Find bond params
        params = []
        key = f'{atoms[bond[0]].atom_type}-{atoms[bond[1]].atom_type}'
        if key in bond_params.keys():
            params = bond_params[key]
        else:
            key_reversed = reverse_param_key(key)
            if key_reversed in bond_params.keys():
                params = bond_params[key_reversed]
            else:
                print('ERROR', key)
        # Write bond params in .itp file
        #         '    2     1     1      0.1529 224262.400'
        file.write(format_val(bond[0]+1,   6   ))
        file.write(format_val(bond[1]+1,   6   ))
        file.write(format_val(params[0], 7   ))
        file.write(format_val(params[1], 7, 4))
        for i in range(2, len(params)):
            file.write(format_val(params[i], 7, 3))
        file.write('\n')
    file.close()
    
def write_angles(atoms, angle_params, filename, skip_long_bonds=False):
    angles = getAllAngles(atoms)
    
    file = open(filename, 'w')
    file.write('[ angles ]\n')
    file.write(';   ai     aj     ak    funct   theta         cth\n')
    for angle in angles:
        if skip_long_bonds:
            delta = max(np.sqrt(sum(np.power((atoms[angle[0]].r - atoms[angle[1]].r), 2))), \
                        np.sqrt(sum(np.power((atoms[angle[1]].r - atoms[angle[2]].r), 2))))
            if delta > 10.0:
                continue
            
        # Find angle params
        params = []
        key = f'{atoms[angle[0]].atom_type}-{atoms[angle[1]].atom_type}-{atoms[angle[2]].atom_type}'
        if key in angle_params.keys():
            params = angle_params[key]
        else:
            key_reversed = reverse_param_key(key)
            if key_reversed in angle_params.keys():
                params = angle_params[key_reversed]
            else:
                # print('ERROR', key)
                continue
        # Write angle params in .itp file
        #         '     1      2      3      1    107.500    502.080'
        file.write(format_val(angle[0]+1,  6   ))
        file.write(format_val(angle[1]+1,  6   ))
        file.write(format_val(angle[2]+1,  6   ))
        file.write(format_val(params[0], 7   ))
        file.write(format_val(params[1], 7, 3))
        for i in range(2, len(params)):
            file.write(format_val(params[i], 7, 3))
        file.write('\n')
    file.close()
    
def write_dihedrals(atoms, dihedral_params, filename, skip_long_bonds=False): 
    dihedrals = getAllDihedrals(atoms)
    
    file = open(filename, 'w')
    file.write('[ dihedrals ]\n')
    file.write(';    i      j      k      l   func   phase     kd     pn\n')
    file.write(';    i      j      k      l   func   c1     c2     c3     c4     c5\n')
    for dihedral in dihedrals:
        if skip_long_bonds:
            delta = max(np.sqrt(sum(np.power((atoms[dihedral[0]].r - atoms[dihedral[1]].r), 2))), \
                        np.sqrt(sum(np.power((atoms[dihedral[1]].r - atoms[dihedral[2]].r), 2))))
            delta = max(delta, \
                        np.sqrt(sum(np.power((atoms[dihedral[2]].r - atoms[dihedral[3]].r), 2))))
            if delta > 10.0:
                continue
            
        # Find dihedral params
        params = []
        key = f'{atoms[dihedral[0]].atom_type}-{atoms[dihedral[1]].atom_type}-{atoms[dihedral[2]].atom_type}-{atoms[dihedral[3]].atom_type}'
        if key in dihedral_params.keys():
            params = dihedral_params[key]
        else:
            key_reversed = reverse_param_key(key)
            if key_reversed in dihedral_params.keys():
                params = dihedral_params[key_reversed]
            else:
                params = [9, 0.00, 0.00, 2]
                
        # Write angle params in .itp file
        file.write(format_val(dihedral[0]+1,   6   ))
        file.write(format_val(dihedral[1]+1,   6   ))
        file.write(format_val(dihedral[2]+1,   6   ))
        file.write(format_val(dihedral[3]+1,   6   ))
        file.write(format_val(params[0],     7   ))
        if params[0] == 9:
            file.write(format_val(params[1], 7, 3))
            file.write(format_val(params[2], 7, 5))
            file.write(format_val(params[3], 7   ))
        else:
            file.write(format_val(params[1], 7, 5))
            for i in range(2, len(params)):
                file.write(format_val(params[i], 5, 5))
        file.write('\n')
    file.close()

def compose_itp_files(itp_filenames, filename):
    file = open(filename, 'w')
    for itp_filename in itp_filenames:
        itp_file = open(itp_filename, 'r')
        for line in itp_file:
            file.write(line)
        file.write('\n')
    file.close()
    
    
    
    
    