import os

adducts_neg = [["[M-H]-", -1.007276466812], ["[M+CH3COO]-", 59.013851]]

names_neg = ["PC 15:0-18:1 d7",
         "LPC 18:1 d7",
         "PE 15:0-18:1 d7",
         "LPE 18:1 d7",
         "PG 15:0-18:1 d7",
         "PI 15:0-18:1 d7",
         "PS 15:0-18:1 d7",
         "SM 18:1-18:1 d9",
         "C15 Ceramide"]

exact_masses_neg = [752.606,
                528.392,
                710.559,
                486.345,
                741.554,
                829.570,
                754.549,
                737.64,
                530.54]

classes_neg = ["PC",
           "LPC",
           "PE",
           "LPE",
           "PG",
           "PI",
           "PS",
           "SM",
           "Cer"]

neg_fragments = [[[241.2486,400],[288.2925,800],[224.069,400],[168.043,400],[152.9958,400]],
                 [[288.2925,800],[224.069,400],[168.043,400],[152.9958,400]],
                 [[241.2486,400],[288.2925,800],[438.2626,50],[196.0375,50],[140.0118,50]],
                 [[288.2925,800],[196.0375,400],[140.0118,400]],
                 [[241.2486,400],[288.2925,800],[516.332,100],[498.322,100],[469.257,100],[451.247,100],[227.033,100],[209.022,100],[152.996,100],[96.97,100],[78.959,100]],
                 [[241.2486,400],[288.2925,800],[604.348,100],[586.338,100],[557.273,100],[539.263,100],[377.2099,100],[315.049,100],[297.038,100],[259.023,100],[241.012,100],[223.001,100],[152.996,100],[96.97,100],[78.959,100]],
                 [[241.2486,400],[288.2925,800],[666.510,100],[442.296,100],[424.285,100], [395.220,100],[377.210,100],[152.996,100],[96.97,100],[78.959,100]],
                 [[722.61,800],[651.56,50]],
                 [[529.53,800]]]

#these retention times are on a different column, so probably wrong..
neg_rts =[6.9,
          4.2,
          7.1,
          4.2,
          5.5,
          5.4,
          5.4,
          6.5,
          7.1]

adducts_pos = [["[M+H]+", 1.007276466812], ["[M+NH4]+", 18.033823]]

names_pos = ["PC 15:0-18:1 d7",
         "LPC 18:1 d7",
         "PE 15:0-18:1 d7",
         "LPE 18:1 d7",
         "TAG 15:0-18:1-15:0 d7",
         "DAG 15:0-18:1 d7",
         "MAG 18:1 d7",
         "CE 18:1 d7",
         "SM 18:1-18:1 d9"]



classes_pos = ["PC",
           "LPC",
           "PE",
           "LPE",
           "TG",
           "DG",
           "MG",
           "CE",
           "SM"]



exact_masses_pos = [752.606,
                528.392,
                710.559,
                486.345,
                811.765,
                587.551,
                363.337,
                657.644,
                737.64]


pos_fragments = [[[184.07,800]],
                 [[184.07,800]],
                 [[299.26,800]],
                 [[299.26,800]],
                 ]

os.chdir("D:\\Johnny\\Lipidomics")

pos_msp = open("pos_lipid_msp.msp", "w")

for i in range(len(names_pos)):
    for j in range(len(adducts_pos)):
        pos_msp.write("NAME: " + names_pos[i] + "\n")
        precursor_mz = round(exact_masses_pos[i] + adducts_pos[j][1], 4)
        pos_msp.write("PRECURSORMZ: " + str(precursor_mz) + "\n")
        pos_msp.write("PRECURSORTYPE: " + adducts_pos[j][0] + "\n")
        pos_msp.write("IONMODE: Positive" + "\n")
        pos_msp.write("COMPOUNDCLASS: " + classes_pos[i] + "\n")
        if i == len(names_neg) - 1 and j == len(adducts_neg) - 1:
            continue
        pos_msp.write("\n")


neg_msp = open("neg_lipid_msp.msp", "w")
for i in range(len(names_neg)):
    for j in range(len(adducts_neg)):
        neg_msp.write("NAME: " + names_neg[i] + "\n")
        precursor_mz = round(exact_masses_neg[i] + adducts_neg[j][1], 4)
        neg_msp.write("PRECURSORMZ: " + str(precursor_mz) + "\n")
        neg_msp.write("PRECURSORTYPE: " + adducts_neg[j][0] + "\n")
        neg_msp.write("RETENTIONTIME: " + str(neg_rts[i]) + "\n")
        neg_msp.write("IONMODE: Negative" + "\n")
        neg_msp.write("COMPOUNDCLASS: " + classes_neg[i] + "\n")
        neg_msp.write("Num Peaks: " + str(len(neg_fragments[i])) + "\n")
        for p in range(len(neg_fragments[i])):
            neg_msp.write(str(neg_fragments[i][p][0]) + " " + str(neg_fragments[i][p][1]) + "\n")
        if i == len(names_neg) - 1 and j == len(adducts_neg) - 1:
            continue
        neg_msp.write("\n")

neg_msp.close()




