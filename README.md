# MD_ZIF7_TEMPO

Aom indices start with 0.

Atomic coordinates are measured in angstroms.

# components

- atom
  - Atom

    name: string, r: np.array(3)
    
    resid_idx: integer, resid: string, atom_idx: integer, atom_type: string

    r_internal: np.zeros(3). Lattice coordinates

    adjacency: list of integers
    
  - mol2ToAtoms(data)
 
    data: output of read_mol2file(), returns list of atoms

  - countAtoms(atoms)
 
    prints number of atoms of each atom_type

# tempo

# zif7

# zif7_tempo

# mdp
