
from sklearn.model_selection import train_test_split
import numpy

import sys

global_NM_fraction_file=sys.argv[1]
output_path=sys.argv[2]
print(global_NM_fraction_file)
print(output_path)
with open(global_NM_fraction_file) as data:
    print(global_NM_fraction_file)
    line = data.readline().strip()
    line = data.readline().strip()
    fraction_dict = {}
    while line:
        peptide = line.split("\t")[1]
        modification = line.split("\t")[2]
        rt = line.split("\t")[4]
        fraction = line.split("\t")[5]
        evalue = line.split("\t")[3]
        fraction = fraction.split(".")[0]        

        if "M" in modification:
            mod_list = modification.split(";")
            for each_mod in mod_list:
                if "M" in each_mod:
                    print("M oxidation")
                    index = int(each_mod.split("@")[1].split("M")[0])
                    peptide = peptide[:index - 1] + "1" + peptide[index:]
        #peptide = "K." + peptide
        #peptide += ".R"

        if fraction in fraction_dict.keys():
            fraction_dict[fraction].append(peptide + "\t" + rt + "\t" + evalue)
        else:
            pep_rt_list = [peptide + "\t" + rt + "\t" + evalue]
            fraction_dict[fraction] = pep_rt_list
        line = data.readline().strip()

    for key in fraction_dict.keys():
        one_file_data = fraction_dict[key]
        data = numpy.array(one_file_data)
        x_train, x_test = train_test_split(data, test_size=0.1)
        train_file = open(
            output_path + key + "_train.txt",
            "w")
        normal_file = open(
            output_path + key + "_normal.txt",
            "w")
        test_file = open(
            output_path + key + "_test.txt",
            "w")
        train_file.write("x\ty\tevalue\n")
        test_file.write("x\ty\tevalue\n")
        normal_file.write("x\ty\tevalue\n")
        for record in x_train:
            train_file.write(record + "\n")
            normal_file.write(record + "\n")
        for record in x_test:
            test_file.write(record + "\n")
            normal_file.write(record + "\n")
