# MD_ZIF7_TEMPO

Atom indices start with 0.

Atom coordinates are measured in angstroms.

# system/components

- atom
  - **Atom**

    name: string, r: np.array(3)
    
    resid_idx: integer, resid: string, atom_idx: integer, atom_type: string

    r_internal: np.zeros(3). Lattice coordinates

    adjacency: list of integers
    
  - copy_atom(atom):

    returns copy of atom
  
  - mol2ToAtoms(data)
 
    data: output of read_mol2file(), returns list of atoms. Usage: atoms = mol2ToAtoms(readMol2File('zif7.mol2'))

  - countAtoms(atoms)
 
    prints number of atoms of each atom_type

  - remove_atoms(atoms, ids_to_remove)

    removes selected atoms from list and corrects adjacency of all atoms

  - remove_non_bonded_atoms
 
    removes atoms with empty adjacency from list and corrects adjacency of all atoms

- read_utils
  
  - read_mol2_file(filename)
 
    returns dictionary where key is section name (AOMS, BONDS, ...) and value is strings of section
    
- write_utils

  - write_gro_file(atoms, filename, a, b, c, alpha=np.pi/2, beta=np.pi/2, gamma=np.pi/2, bounds_a=[0,1], bounds_b=[0,1], bounds_c=[0,1], scale=1.0)
 
  - write_mol2_file(atoms, filename, a, b, c, alpha, beta, gamma, skip_long_bonds=False)
 
  - write_atoms(atoms, charge, mass, resid, filename)
 
  - write_bonds(atoms, bond_params, filename)
 
  - write_angles(atoms, angle_params, filename)
 
  - write_dihedrals(atoms, dihedral_params, filename)
 
  - compose_itp_files(itp_filenames, filename)

- get_topology_params

- pbc

  - calculate_lattice_vectors(a, b, c, alpha, beta, gamma)
 
  - compare_periodic(r1, r2, p_a, p_b, p_c)
 
  - apply_pbc(atoms, a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
# system/tempo

# system/zif7

# system/zif7_tempo

# mdp
