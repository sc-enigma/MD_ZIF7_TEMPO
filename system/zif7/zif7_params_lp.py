def get_zif7_params_lp():
    # [atom_type] = m
    mass = {}
    mass['intraFF_Zn'] = 65.360
    mass['intraFF_N']  = 14.007
    mass['intraFF_C1'] = 12.011
    mass['intraFF_C2'] = 12.011
    mass['intraFF_C5'] = 12.011
    mass['intraFF_C6'] = 12.011
    mass['intraFF_H1'] =  1.008
    mass['intraFF_H5'] =  1.008
    mass['intraFF_H6'] =  1.008
    mass['repl_O'] =     16.000
    mass['repl_H'] =      1.008
 
    # [atom_type] = q
    charge = {}
    charge['intraFF_Zn'] = 0.7828
    charge['intraFF_N']  =-0.4502
    charge['intraFF_C1'] = 0.2432
    charge['intraFF_C2'] = 0.1722
    charge['intraFF_C5'] =-0.0881
    charge['intraFF_C6'] =-0.1982
    charge['intraFF_H1'] = 0.0988
    charge['intraFF_H5'] = 0.0988
    charge['intraFF_H6'] = 0.0988
    charge['repl_O']     =-0.5893
    charge['repl_H']     = 0.4055

    # [atom_type-atom_type] = [funct, r0, k]
    bond_params = {}
    bond_params['intraFF_C1-intraFF_H1'] = [1, 0.1090, 154544.408]
    bond_params['intraFF_C1-intraFF_N']  = [1, 0.1362, 146958.816]
    bond_params['intraFF_C2-intraFF_N']  = [1, 0.1399, 114733.648]
    bond_params['intraFF_C2-intraFF_C2'] = [1, 0.1447, 120637.272]
    bond_params['intraFF_C2-intraFF_C5'] = [1, 0.1357, 146561.336]
    bond_params['intraFF_C5-intraFF_H5'] = [1, 0.1090, 151360.384]
    bond_params['intraFF_C5-intraFF_C6'] = [1, 0.1408, 153523.512]
    bond_params['intraFF_C6-intraFF_H6'] = [1, 0.1090, 150962.904]
    bond_params['intraFF_C6-intraFF_C6'] = [1, 0.1408, 142247.632]
    bond_params['intraFF_Zn-intraFF_N']  = [3, 0.1090, 109.0350, 19.8000]
    bond_params['repl_O-intraFF_Zn']     = [1, 0.1391, 476976.000]
    bond_params['repl_O-repl_H']         = [1, 0.1933, 462750.400]

    # [atom_type-atom_type-atom_type] = [funct, angle, k]
    angle_params = {}
    angle_params['intraFF_N-intraFF_C1-intraFF_H1']  = [1, 122.77, 227.610]
    angle_params['intraFF_N-intraFF_C1-intraFF_N']   = [1, 111.90, 473.880]
    angle_params['intraFF_N-intraFF_C2-intraFF_C5']  = [1, 130.92, 628.939]
    angle_params['intraFF_N-intraFF_C2-intraFF_C2']  = [1, 106.62, 510.406]
    angle_params['intraFF_C1-intraFF_N-intraFF_C2']  = [1, 106.39, 598.479]
    angle_params['intraFF_C2-intraFF_C5-intraFF_H5'] = [1, 121.52, 251.542]
    angle_params['intraFF_C2-intraFF_C5-intraFF_C6'] = [1, 115.95, 658.436]
    angle_params['intraFF_C5-intraFF_C6-intraFF_C6'] = [1, 120.50, 796.968]
    angle_params['intraFF_C5-intraFF_C6-intraFF_H6'] = [1, 119.17, 251.458]
    angle_params['intraFF_C6-intraFF_C6-intraFF_H6'] = [1, 119.15, 253.341]
    angle_params['intraFF_C6-intraFF_C5-intraFF_H5'] = [1, 121.48, 267.274]
    angle_params['intraFF_C2-intraFF_C2-intraFF_C5'] = [1, 119.82, 610.027]
    angle_params['intraFF_N-intraFF_Zn-intraFF_N']   = [1, 109.42, 134.055]
    angle_params['intraFF_Zn-intraFF_N-intraFF_C1']  = [1, 126.76, 206.020]
    angle_params['intraFF_Zn-intraFF_N-intraFF_C2']  = [1, 127.94, 214.095]
    angle_params['intraFF_N-intraFF_Zn-repl_O']      = [1, 110.97, 76.3000]
    angle_params['intraFF_Zn-repl_O-repl_H']         = [1, 112.68, 76.3000]
    angle_params['repl_H-repl_O-repl_H']             = [1, 109.47, 383, 109.47, 383]

    # [atom_type-atom_type-atom_type-atom_type] = [funct, angle, k, n]        - periodic
    # [atom_type-atom_type-atom_type-atom_type] = [funct, c1, c2, c3, c4, c5] - fourier
    dihedral_params = {}
    dihedral_params['intraFF_N-intraFF_C2-intraFF_C5-intraFF_C6']  = [9, 180.00, 3.09616, 2]
    dihedral_params['intraFF_H1-intraFF_C1-intraFF_N-intraFF_C2']  = [9, 180.00, 15.31344, 2]
    dihedral_params['intraFF_C1-intraFF_N-intraFF_C2-intraFF_C2']  = [9, 180.00, 5.27184, 2]
    dihedral_params['intraFF_H5-intraFF_C5-intraFF_C2-intraFF_N']  = [9, 180.00, 1.38072, 2]
    dihedral_params['intraFF_C2-intraFF_C2-intraFF_C5-intraFF_H5'] = [9, 180.00, 1.88280, 2]
    dihedral_params['intraFF_C5-intraFF_C2-intraFF_C2-intraFF_C5'] = [9, 180.00, 3.93296, 2]
    dihedral_params['intraFF_H5-intraFF_C5-intraFF_C6-intraFF_C6'] = [9, 180.00, 1.79912, 2]
    dihedral_params['intraFF_H6-intraFF_C6-intraFF_C6-intraFF_C5'] = [9, 180.00, 3.59824, 2]
    dihedral_params['intraFF_C5-intraFF_C6-intraFF_C6-intraFF_C5'] = [9, 180.00, 9.03744, 2]
    dihedral_params['intraFF_N-intraFF_C2-intraFF_C2-intraFF_N']   = [9, 180.00, 9.37216, 2]
    dihedral_params['intraFF_N-intraFF_C2-intraFF_C2-intraFF_C5']  = [9, 180.00, 3.76560, 2]
    dihedral_params['intraFF_C1-intraFF_N-intraFF_C2-intraFF_C5']  = [9, 180.00, 4.39320, 2]
    dihedral_params['intraFF_C2-intraFF_C5-intraFF_C6-intraFF_H6'] = [9, 180.00, 10.75288, 2]
    dihedral_params['intraFF_C6-intraFF_C5-intraFF_C2-intraFF_C2'] = [9, 180.00, 2.25936, 2]
    dihedral_params['intraFF_C6-intraFF_C6-intraFF_C5-intraFF_C2'] = [9, 180.00, 4.26768, 2]
    dihedral_params['intraFF_N-intraFF_C1-intraFF_N-intraFF_C2']   = [9, 180.00, 14.85320, 2]
    dihedral_params['intraFF_C1-intraFF_N-intraFF_Zn-intraFF_N']   = [5, 18.82800, 15.39712, 0.54392, 0.00000]
    dihedral_params['intraFF_C2-intraFF_N-intraFF_Zn-intraFF_N']   = [5, 17.27992, 11.33864, 0.16736, 0.00000]

















    return mass, charge, bond_params, angle_params, dihedral_params