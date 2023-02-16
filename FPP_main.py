import os, sys, numpy, tabulate
from Bio.SVDSuperimposer import SVDSuperimposer
from FPP_fun import compute_pro_parameters

PATH_1 = os.getcwd() + "\x5c" + sys.argv[1]
PATH_2 = os.getcwd() + "\x5c" + sys.argv[2]
# PATH_1 = "C:/Users/aless/PycharmProjects/pdbparser/3epo.pdb"
# PATH_2 = "C:/Users/aless/PycharmProjects/pdbparser/1elo.pdb"

#   Script can be used by running the code with files as arguments
#   or by plugging the absolute path of both files in the code.


pro_mass_1, pro_volume_1, center_mass_coord_1, distance_FL_1, atom_pos_1,\
    counter_atom_errors_1, counter_residue_errors_1\
    = compute_pro_parameters(PATH_1)
#   0:pro_mass, 1:pro_volume, 2:center_mass_coord, 3:distance_FL, 4:atom_pos
#   5:counter_atom_errors, 6:counter_residue_errors
pro_mass_2, pro_volume_2, center_mass_coord_2, distance_FL_2, atom_pos_2,\
    counter_atom_errors_2, counter_residue_errors_2\
    = compute_pro_parameters(PATH_2)


if counter_atom_errors_1 > 0 or counter_residue_errors_1 > 0:
    print("WARNING: check file", sys.argv[1], "there might be some errors")
if counter_atom_errors_2 > 0 or counter_residue_errors_2 > 0:
    print("WARNING: check file", sys.argv[2], "there might be some errors")


print("\n")

print(tabulate.tabulate([
                ["Protein Mass", pro_mass_1, pro_mass_2],
                ["Protein Volume", pro_volume_1, pro_volume_2],
                ["Center of Mass", center_mass_coord_1, center_mass_coord_2],
                ["End-to-End distance", distance_FL_1, distance_FL_2],
                ["Atoms Not Counted", counter_atom_errors_1, counter_atom_errors_2],
                ["Residues Cot Counted", counter_residue_errors_1, counter_residue_errors_2]
                ],
               headers=[sys.argv[1], sys.argv[2]]
                ))

atom_pos_1 = numpy.array(atom_pos_1, dtype='f')
atom_pos_2 = numpy.array(atom_pos_2, dtype='f')

print("\n")

try:
    superimposed = SVDSuperimposer()
    superimposed.set(atom_pos_1, atom_pos_2)
    superimposed.run()
    rmsd = superimposed.get_rms()
    print("The computed RMSD is: ", round(rmsd, 6), "Å")
except:
    print("Couldn't compute RMSD; proteins have different backbone length")

# from FPP_RMSD import compute_RMSD
#
# rmsd = compute_RMSD(atom_pos_1, atom_pos_2)
# if rmsd is not None:
#     print("The computed RMSD is: ", round(rmsd, 6), "Å")
# else:
#     print("Couldn't compute RMSD; proteins have different backbone length")
