def compute_pro_parameters(PATH):
    counter_atom_errors = int()
    counter_residue_errors = int()

    residue_mass_dict = {
        'ALA': '71.03711',
        'CYS': '103.00919',
        'ASP': '115.02694',
        'GLU': '129.04259',
        'PHE': '147.06841',
        'GLY': '57.02146',
        'HIS': '137.05891',
        'ILE': '113.08406',
        'LYS': '128.09496',
        'LEU': '113.08406',
        'MET': '131.04049',
        'ASN': '114.04293',
        'PRO': '97.05276',
        'GLN': '128.05858',
        'ARG': '156.10111',
        'SER': '87.03203',
        'THR': '101.04768',
        'VAL': '99.06841',
        'TRP': '186.07931',
        'TYR': '163.06333'
    }
#   Hydrogen is always considered for computing the protein mass
    pro_mass = float()
    AA = None

    atom_mass_dict = {
        'C': '12.011',
        'O': '15.9994',
        'N': '14.0067',
        'S': '32.06',
        'H': '1.00784'
    }
#   Hydrogen are/aren't considered for computing the center of mass depending on the PDB file:
#   if hydrogen are present they will be considered, else they will not
    sum_of_mass_x = float()
    sum_of_mass_y = float()
    sum_of_mass_z = float()
    atom_mass_sum = float()

    flag = True

    atom_pos = list()

    pro_density = 0.8129891
#   Density expressed as Dalton/cubic Angstrom (converted from average protein density of 1.35 g/cm³)

    with open(PATH) as pdb:
        for line in pdb.readlines():
            if line.startswith("ATOM"):
                x = line.split()

                if len(x) < 12:
                    counter_atom_errors += 1
                    continue

                if x[3] != AA:
                    if residue_mass_dict.get(x[3]) is not None:
                        AA = x[3]
                        pro_mass += float(residue_mass_dict[x[3]])
                    else:
                        counter_residue_errors += 1

                atom_mass_sum += float(atom_mass_dict[x[11]])
                sum_of_mass_x += float(x[6]) * float(atom_mass_dict[x[11]])
                sum_of_mass_y += float(x[7]) * float(atom_mass_dict[x[11]])
                sum_of_mass_z += float(x[8]) * float(atom_mass_dict[x[11]])

                if flag:
                    first_x = float(x[6])
                    first_y = float(x[7])
                    first_z = float(x[8])
                    flag = False
                first_last_x = float(x[6]) - first_x
                first_last_y = float(x[7]) - first_y
                first_last_z = float(x[8]) - first_z

                if x[2] == "CA":
                    atom_pos.append(list([float(x[6]), float(x[7]), float(x[8])]))

    pro_volume = pro_mass / pro_density
    if pro_volume > 1000:
        pro_volume = pro_volume / 1000
        pro_volume = str(round(pro_volume, 3)) + " nm³"
    else:
        pro_volume = str(round(pro_volume, 3)) + " Å³"

    if pro_mass > 1000:
        pro_mass = pro_mass / 1000
        pro_mass = str(round(pro_mass, 3)) + " KDa"
    else:
        pro_mass = str(round(pro_mass, 3)) + " Da"

    center_mass_x = round(sum_of_mass_x / atom_mass_sum, 3)
    center_mass_y = round(sum_of_mass_y / atom_mass_sum, 3)
    center_mass_z = round(sum_of_mass_z / atom_mass_sum, 3)
    center_mass_coord = (center_mass_x, center_mass_y, center_mass_z)

    distance_FL = round(((first_last_x ** 2) + (first_last_y ** 2) + (first_last_z ** 2)) ** (1 / 2), 3)

    return pro_mass, pro_volume, center_mass_coord, distance_FL, atom_pos, counter_atom_errors, counter_residue_errors
