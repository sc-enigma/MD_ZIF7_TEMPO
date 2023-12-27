# MD_ZIF7_TEMPO

Atom indices start with 0.

Atom coordinates are measured in angstroms.

# system/components

- **atom.py**
  - **Atom**

      name: string, r: np.array(3)
    
      resid_idx: integer, resid: string, atom_idx: integer, atom_type: string

      r_internal: np.zeros(3). Lattice coordinates

      adjacency: list of integers
    
  - **copy_atom(atom)**

    returns copy of atom
  
  - **mol2ToAtoms(data)**
 
    data: output of read_mol2file(), returns list of atoms. Usage: atoms = mol2ToAtoms(readMol2File('zif7.mol2'))

  - **countAtoms(atoms)**
 
    prints number of atoms of each atom_type

  - **remove_atoms(atoms, ids_to_remove)**

    removes selected atoms from list and corrects adjacency of all atoms

  - **remove_non_bonded_atoms**
 
    removes atoms with empty adjacency from list and corrects adjacency of all atoms

- **read_utils.py**
  
  - read_mol2_file(filename)
 
    returns dictionary where key is section name (AOMS, BONDS, ...) and value is strings of section
    
- **write_utils.py**

  - **write_gro_file(atoms, filename, a, b, c, alpha=np.pi/2, beta=np.pi/2, gamma=np.pi/2, bounds_a=[0,1], bounds_b=[0,1], bounds_c=[0,1], scale=1.0)**
 
  - **write_mol2_file(atoms, filename, a, b, c, alpha, beta, gamma, skip_long_bonds=False)**
 
  - **write_atoms(atoms, charge, mass, resid, filename)**
    Write GROMACS itp file. Charge and mass are dictionaries with key = atomtype and value = double
 
  - **write_bonds(atoms, bond_params, filename)**
    Write GROMACS itp file. bond_params is dictionaries with key = atomtype and value = list of bond params
 
  - **write_angles(atoms, angle_params, filename)**
    Write GROMACS itp file. angle_params is dictionaries with key= atomtype and value = list of angle params
 
  - **write_dihedrals(atoms, dihedral_params, filename)**
    Write GROMACS itp file. dihedral_params is dictionaries with key = atomtype and value = list of dihedral params
 
  - **compose_itp_files(itp_filenames, filename)**
    Unites several files in one. itp_filenames is a list of input filenames

- **get_topology_params.py**
    - **getAllBonds(atoms)**
      Returns list of all unique bonds in atoms. Bond is a list consists of 2 atomic indices
 
    - **getAllAngles(atoms)**
      Returns list of all unique angles in atoms. Angle is a list consists of 3 atomic indices
 
    - **getAllDihedrals(atoms)**
      Returns list of all unique dihedrals in atoms. Dihedral is a list consists of 4 atomic indices

- **pbc.py**

  - **calculate_lattice_vectors(a, b, c, alpha, beta, gamma)**
    Calculates lattice vectors
 
  - **compare_periodic(r1, r2, p_a, p_b, p_c)**
    Compares 2 lattice vectors. p_a, p_b, p_c are periods (default 0)
 
  - **apply_pbc(atoms, a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)**
    Calculates lattice vectors for each atom. All atoms for which at least one lattice vector coordinate is not in the corresponding interval are removed: bounds_a, bounds_b, bounds_c. For the remaining atoms, the adjacency is adjusted to ensure correct periodic boundary conditions.
    
# system/tempo
Topology of TEMPO

# system/zif7
Topology of ZIF7 + scripts for topology preparation

# system/zif7_tempo
Topology of TEMPO in ZIF7 + scripts for topology preparation

# system/zif7_tempo_removed
Topology of TEMPO in ZIF7 with removed linkers + scripts for topology preparation

# system/zif7_tempo_CO2
Topology of TEMPO in ZIF7 with absorbed CO2 molecules + scripts for topology preparation

# mdp

- em.mdp
- eql.mdp
- prod.mdp
