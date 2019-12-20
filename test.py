# import glob
#
# read_files = glob.glob("sequences/clustalo/only_sequences/sequences/*.fasta")
#
# with open("sequences/clustalo/only_sequences/all_sequences.fasta", "w") as outfile:
#     for f in read_files:
#         with open(f, "r") as infile:
#             outfile.write(infile.read() + "\n" + "\n")

# from shutil import copyfile, copy
# import glob
#
# read_files_minus = glob.glob("sequences/ssRNA_minus/*fasta")
#
# listeminus = []
#
# for file in read_files_minus:
#     with open(file, "r") as f:
#         if "complete cds" in f.readline():
#             listeminus.append(file)
#
# for i in listeminus:
#     copy(i, "sequences/clustalo/only_cds/sequences/")
#
# read_files_plus = glob.glob("sequences/ssRNA_plus/*fasta")
#
# listeplus = []
#
# for file in read_files_plus:
#     with open(file, "r") as f:
#         if "complete cds" in f.readline():
#             listeplus.append(file)
#
# for i in listeplus:
#     copy(i, "sequences/clustalo/only_cds/sequences/")

#np.savetxt("test1.csv", names, delimiter=",")

#array = array[:20, :20]
#names = names[:20]

# im = plt.imshow(array, vmin=0, vmax=1)
# ax=plt.gca()
#
# ax.set_xticks(np.arange(len(names)))
# ax.set_yticks(np.arange(len(names)))
#
# ax.set_xticklabels(names)
# ax.set_yticklabels(names)
#
# plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
#
# plt.colorbar()
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
#
# with open("sequences/clustalo/only_cds/distmat_cds.mat", "r") as file:
#     size = int(file.readline())
#     array = np.zeros((size, size))
#     names = ['']*size
#
#     for i, line in enumerate(file):
#         names[i] = line[0:11].strip()
#         stripped = line[12:].strip()
#         floats = [float(x) for x in stripped.split()]
#         floats_to_array = np.array(floats)
#         array[i] = floats_to_array
#
# np.savetxt("test.csv", array, delimiter=",")

import numpy as np
import csv
import matplotlib.pyplot as plt

with open("sequences/clustalo/all_sequences/distmat_all.mat", "r") as file:
    size = int(file.readline())
    array = np.zeros((size, size))
    names = ['']*size

    for i, line in enumerate(file):
        names[i] = line[0:11].strip()
        stripped = line[12:].strip()
        floats = [float(x) for x in stripped.split()]
        floats_to_array = np.array(floats)
        array[i] = floats_to_array

np.savetxt("sequences/clustalo/all_sequences/test.csv", array, delimiter=",")

with open("sequences/clustalo/all_sequences/test2.csv", 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(names)

with open("sequences/clustalo/all_sequences/test1.csv", "w") as outfile:
    with open("sequences/clustalo/all_sequences/test.csv", "r") as infile:
        i = 0
        for line in infile:
            outfile.write("\"" + names[i] + "\"" + "," + line)
            i = i + 1