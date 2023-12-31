import json

GERDA_BEGe_phII = {
    "isotope": ["42K", "40K", "60Co", "60Co", "228Ac", "228Ac", "228Ac", "228Ac", "212Bi", "212Pb", "208Tl", "208Tl", "208Tl", "214Pb", "214Pb", "214Pb", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "e+e-", "85Kr", "65Zn"],
    "energy": [1525, 1461, 1333, 1173, 911, 969, 965, 338, 727, 239, 2615, 583, 861, 295, 242, 352, 609, 1765, 1238, 2204, 1378, 2448, 1120, 511, 514, 1125],
    "raw": {
        "mode": [76.6, 49.3, 0.8, 1.4, 3.0, 1.7, 0.6, 2.2, 1.5, 13.2, 0.9, 2.8, 0.7, 7.3, 3.9, 10.8, 4.9, 1.3, 1.3, 0.5, 0.6, 0.2, 1.1, 2.6, 5.7, 0.6],
        "low_68": [1.6, 1.4, 0.4, 0.5, 0.7, 0.7, 0, 1.2, 0.7, 2.8, 0.1, 0.5, 0.6, 1.6, 0, 1.2, 0.7, 0.3, 0.5, 0.1, 0, 0.1, 0.5, 0.7, 0.8, 0],
        "upp_68": [1.6, 1.2, 0.5, 0.7, 0.7, 0.6, 0, 1.2, 0.7, 2.1, 0.2, 1.0, 0.4, 1.6, 0, 1.5, 0.9, 0.2, 0.5, 0.1, 0, 0.1, 0.6, 0.8, 0.9, 0]
    },
    "LArAC": {
        "mode": [16.1, 47.7, 0.6, 0.9, 1.4, 0.9, 0.3, 1.0, 1.0, 10.5, 0.1, 0.7, 1.0, 5.1, 2.4, 7.5, 1.2, 0.6, 0.7, 0.4, 0.4, 0.09, 0.9, 0.4, 3.6, 0.4],
        "low_68": [0.8, 1.2, 0, 0.5, 0.5, 0.5, 0, 0.9, 0.6, 1.6, 0, 0.5, 0, 0.9, 1.4, 0.8, 0, 0.2, 0.4, 0.1, 0, 0.05, 0, 0, 0.7, 0],
        "upp_68": [0.7, 1.3, 0, 0.4, 0.7, 0.6, 0, 0.6, 0.6, 1.5, 0, 0.4, 0, 1.4, 1.2, 1.2, 0, 0.1, 0.4, 0.1, 0, 0.07, 0, 0, 0.6, 0]
    },
    "LArC": {
        "mode": [60.3, 1.5, 0.7, 0.7, 1.5, 0.7, 0.4, 1.3, 0.5, 2.8, 0.8, 2.4, 0.6, 2.0, 1.9, 3.4, 4.5, 0.7, 0.7, 0.13, 0.5, 0.1, 0.9, 3.7, 2.0, 0.5],
        "low_68": [1.1, 0.3, 0.3, 0.3, 0.4, 0.3, 0.3, 0.9, 0.3, 1.8, 0.2, 0.5, 0.4, 1.1, 0, 0.9, 0.7, 0.2, 0.4, 0.06, 0, 0.05, 0.3, 0.6, 0.4, 0],
        "upp_68": [1.7, 0.3, 0.3, 0.3, 0.4, 0.4, 0.3, 0.7, 0.4, 1.5, 0.1, 0.5, 0.3, 1.3, 0, 0.7, 0.3, 0.2, 0.3, 0.1, 0, 0.08, 0.4, 0.6, 0.6, 0]
    }
}



GERDA_BEGe_phII_plus = {
    "isotope": ["42K", "40K", "60Co", "60Co", "228Ac", "228Ac", "228Ac", "228Ac", "212Bi", "212Pb", "208Tl", "208Tl", "208Tl", "214Pb", "214Pb", "214Pb", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "e+e-", "85Kr", "65Zn"],
    "energy": [1525, 1461, 1333, 1173, 911, 969, 965, 338, 727, 239, 2615, 583, 861, 295, 242, 352, 609, 1765, 1238, 2204, 1378, 2448, 1120, 511, 514, 1125],
    "raw": {
        "mode": [75.5, 44.5, 1.7, 0.9, 2.8, 0.6, 1.6, 2.9, 1.9, 10.4, 0.8, 1.0, 1.0, 7.6, 5.3, 7.3, 4.9, 1.0, 0.6, 0.3, 1.0, 0.3, 0.6, 3.1, 4.9, 0.6],
        "low_68": [1.6, 1.6, 0.5, 0.6, 0.7, 0.5, 0.7, 1.3, 0, 3.0, 0.2, 0.8, 0, 1.9, 0, 1.3, 0.8, 0.2, 0.5, 0.1, 0.4, 0.4, 0.5, 0.9, 0.9, 0.4],
        "upp_68": [2.2, 1.4, 0.6, 0.5, 0.9, 0.7, 0.7, 1.3, 0, 2.7, 0.2, 0.6, 0, 1.8, 0, 1.3, 1.1, 0.2, 0.4, 0.1, 0.5, 0.5, 0.6, 0.8, 1.0, 0.6]
    },
    "LArAC": {
        "mode": [14.7, 43.1, 1.1, 0.5, 1.3, 0.8, 1.2, 2.7, 1.3, 5.1, 0.1, 0.5, 0.5, 3.7, 1.6, 4.7, 1.1, 0.6, 0.6, 0.2, 0.7, 0.24, 0.7, 1.2, 2.1, 0.8],
        "low_68": [0.8, 1.3, 0.4, 0, 0.6, 0, 0.5, 0.8, 0, 1.3, 0, 0, 0, 1.1, 0, 0.9, 0.6, 0.2, 0, 0.1, 0.3, 0.08, 0, 0, 0.7, 0.4],
        "upp_68": [0.9, 1.6, 0.4, 0, 0.6, 0, 0.6, 1.2, 0, 1.7, 0, 0, 0, 1.3, 0, 0.9, 0.7, 0.2, 0, 0.1, 0.4, 0.14, 0, 0, 0.6, 0.6]
    },
    "LArC": {
        "mode": [61.1, 1.2, 0.7, 1.1, 1.8, 1.1, 0.7, 1.7, 0.6, 3.3, 0.8, 1.5, 0.3, 3.8, 4.2, 2.5, 3.8, 0.5, 0.9, 0.09, 0.4, 0.05, 0.9, 3.1, 3.0, 0.5],
        "low_68": [1.8, 0.3, 0.4, 0.4, 0.6, 0.3, 0.3, 0, 0.4, 2.3, 0.2, 0.5, 0.2, 1.3, 0, 0.8, 0.6, 0.2, 0, 0.1, 0.3, 0.07, 0.4, 0.5, 0.6, 0],
        "upp_68": [1.6, 0.3, 0.3, 0.5, 0.4, 0.5, 0.4, 0, 0.4, 1.5, 0.2, 0.5, 0.4, 1.3, 0, 0.9, 0.6, 0.2, 0, 0.06, 0.2, 0.05, 0.4, 0.7, 0.6, 0]
    }
}

GERDA_COAX_phII = {
    "isotope": ["42K", "40K", "60Co", "60Co", "228Ac", "228Ac", "228Ac", "228Ac", "212Bi", "212Pb", "208Tl", "208Tl", "208Tl", "214Pb", "214Pb", "214Pb", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "e+e-", "85Kr", "65Zn"],
    "energy": [1525, 1461, 1333, 1173, 911, 969, 965, 338, 727, 239, 2615, 583, 861, 295, 242, 352, 609, 1765, 1238, 2204, 1378, 2448, 1120, 511, 514, 1125],
    "raw": {
        "mode": [109.7, 59.5, 2.0, 2.2, 0.7, 1.1, 0.5, 2.2, 0.8, 3.8, 1.1, 2.2, 1.4, 3.5, 5.6, 6.2, 4.6, 1.8, 1.2, 0.4, 0.8, 0.2, 2.3, 5.2, 5.7, 1.6],
        "low_68": [2.3, 1.8, 0.6, 0.8, 0.5, 0.7, 0.5, 1.7, 0.6, 2.8, 0.2, 0.6, 0, 2.3, 0, 1.7, 0.9, 0.3, 0, 0.1, 0, 0.1, 0.6, 1.1, 1.4, 0.7],
        "upp_68": [1.7, 1.4, 0.6, 0.6, 0.7, 0.6, 0.6, 1.2, 0.7, 2.4, 0.3, 1.0, 0, 2.1, 0, 1.8, 1.0, 0.2, 0, 0.1, 0, 0.1, 0.7, 1.1, 0.8, 0.6]
    },
    "LArAC": {
        "mode": [23.5, 58.5, 0.3, 1.1, 1.0, 0.4, 0.8, 1.7, 1.4, 2.3, 0.2, 1.6, 1.0, 3.1, 1.4, 4.7, 0.9, 0.9, 0.6, 0.19, 0.3, 0.15, 1.3, 0.8, 4.4, 0.7],
        "low_68": [1.1, 1.5, 0.2, 0, 0, 0.3, 0.6, 1.2, 0, 1.5, 0, 0, 0, 1.8, 0, 1.1, 0.6, 0.2, 0, 0.09, 0, 0.07, 0.5, 0.6, 0.7, 0.4],
        "upp_68": [0.9, 1.6, 0.4, 0, 0, 0.5, 0.5, 1.5, 0, 2.1, 0, 0, 0, 1.8, 0, 1.8, 0.6, 0.2, 0, 0.10, 0, 0.08, 0.6, 0.7, 0.9, 0.6]
    },
    "LArC": {
        "mode": [85.4, 0.9, 1.8, 1.7, 0.9, 0.5, 0.5, 1.6, 0.5, 0.7, 1.2, 1.7, 0.5, 2.5, 0.9, 0.9, 4.0, 0.9, 0.7, 0.15, 0.5, 0.1, 1.0, 4.8, 0.4, 0.8],
        "low_68": [1.3, 0.3, 0.4, 0.4, 0.3, 0.3, 0, 0, 0.3, 0, 0.2, 0.4, 0.4, 0, 0, 0.7, 0.7, 0.2, 0.5, 0.07, 0.3, 0, 0.4, 0.8, 0.4, 0.4],
        "upp_68": [2.2, 0.3, 0.4, 0.5, 0.5, 0.4, 0, 0, 0.4, 0, 0.2, 0.5, 0.3, 0, 0, 0.8, 0.5, 0.2, 0.3, 0.11, 0.3, 0, 0.4, 0.6, 0.5, 0.3]
    }
}



GERDA_COAX_phII_plus = {
    "isotope": ["42K", "40K", "60Co", "60Co", "228Ac", "228Ac", "228Ac", "228Ac", "212Bi", "212Pb", "208Tl", "208Tl", "208Tl", "214Pb", "214Pb", "214Pb", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "e+e-", "85Kr", "65Zn"],
    "energy": [1525, 1461, 1333, 1173, 911, 969, 965, 338, 727, 239, 2615, 583, 861, 295, 242, 352, 609, 1765, 1238, 2204, 1378, 2448, 1120, 511, 514, 1125],
    "raw": {
        "mode": [96.4, 50.0, 2.1, 2.1, 0.9, 2.1, 1.8, 4.3, 2.7, 4.3, 1.1, 4.3, 1.4, 6.2, 8.3, 8.9, 2.7, 1.8, 1.5, 0.3, 0.9, 0.15, 1.3, 0.8, 4.4, 0.7],
        "low_68": [2.7, 1.9, 0.8, 0, 0.7, 0, 0, 0, 1.2, 3.9, 0.3, 1.2, 0, 2.8, 0, 2.8, 1.2, 0.4, 0, 0.1, 0, 0, 0.9, 1.6, 1.4, 0],
        "upp_68": [2.8, 2.5, 0.9, 0, 1.2, 0, 0, 0, 1.3, 2.9, 0.3, 1.5, 0, 3.7, 0, 2.6, 1.2, 0.4, 0, 0.3, 0, 0, 0.7, 1.3, 1.6, 0]
    },
    "LArAC": {
        "mode": [16.3, 49.7, 1.7, 1.6, 2.0, 1.4, 1.2, 2.0, 1.7, 4.9, 0.3, 0.9, 1.1, 5.2, 6.8, 7.2, 0.9, 0.7, 1.5, 0.09, 1.0, 0.1, 0.8, 0.9, 2.2, 1.0],
        "low_68": [1.3, 1.8, 1.2, 0, 0, 0, 0, 1.5, 1.2, 2.9, 0, 0.8, 0, 2.2, 0, 2.6, 0, 0.3, 0, 0.07, 0, 0, 0, 0.8, 1.0, 0],
        "upp_68": [1.4, 2.4, 0.8, 0, 0, 0, 0, 2.0, 0.8, 3.0, 0, 1.0, 0, 2.9, 0, 1.9, 0, 0.3, 0, 0.14, 0, 0, 0, 0.8, 1.3, 0]
    },
    "LArC": {
        "mode": [78.5, 0.9, 1.0, 1.3, 1.1, 0.9, 0.4, 1.1, 1.1, 5.3, 0.9, 3.6, 1.0, 4.1, 5.0, 2.4, 3.6, 1.0, 0.6, 0.25, 0.5, 0.2, 1.7, 3.0, 2.2, 1.1],
        "low_68": [2.1, 0.4, 0.5, 0, 0.6, 0.6, 0.3, 0, 0.6, 0, 0.2, 0.7, 0, 0, 0, 1.4, 0.6, 0.3, 0, 0.14, 0, 0, 0.7, 0.8, 0.8, 0],
        "upp_68": [2.9, 0.4, 0.6, 0, 0.5, 0.4, 0.6, 0, 0.6, 0, 0.4, 0.9, 0, 0, 0, 1.2, 0.9, 0.3, 0, 0.17, 0, 0, 0.5, 1.1, 1.0, 0]
    }
}

GERDA_IC_phII_plus = {
    "isotope": ["42K", "40K", "60Co", "60Co", "228Ac", "228Ac", "228Ac", "228Ac", "212Bi", "212Pb", "208Tl", "208Tl", "208Tl", "214Pb", "214Pb", "214Pb", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "e+e-", "85Kr", "65Zn"],
    "energy": [1525, 1461, 1333, 1173, 911, 969, 965, 338, 727, 239, 2615, 583, 861, 295, 242, 352, 609, 1765, 1238, 2204, 1378, 2448, 1120, 511, 514, 1125],
    "raw": {
        "mode": [120.9, 99.9, 2.0, 1.2, 2.9, 2.9, 2.3, 2.6, 1.5, 13.8, 1.8, 1.6, 2.6, 5.4, 8.0, 8.7, 4.9, 2.5, 0.9, 1.2, 1.2, 0.2, 3.8, 6.1, 5.1, 10.5],
        "low_68": [3.3, 2.7, 0, 0, 1.4, 1.2, 0, 1.7, 1.1, 3.6, 0.4, 1.1, 0, 1.8, 4.0, 1.7, 1.1, 0.6, 0.8, 0.3, 0.8, 0, 1.3, 1.5, 1.5, 1.5],
        "upp_68": [4.4, 4.3, 0, 0, 1.1, 1.2, 0, 2.1, 1.1, 4.1, 0.5, 1.2, 0, 3.5, 3.8, 2.3, 1.5, 0.5, 0.7, 0.5, 0.8, 0, 1.4, 1.8, 1.4, 1.6]
    },
    "LArAC": {
        "mode": [19.1, 97.9, 1.1, 0.9, 2.3, 1.8, 1.6, 4.3, 1.0, 12.0, 0.1, 1.9, 2.1, 5.8, 3.9, 6.9, 1.4, 0.8, 0.9, 0.3, 0.6, 0.2, 1.9, 1.6, 2.0, 9.9],
        "low_68": [1.3, 3.1, 0, 0, 0, 1.0, 0, 1.6, 0.8, 2.7, 0.1, 0, 0, 1.7, 0, 1.5, 0.8, 0.3, 0.6, 0.2, 0.5, 0, 1.0, 1.0, 0.9, 1.4],
        "upp_68": [1.8, 3.7, 0, 0, 0, 1.0, 0, 1.6, 1.0, 2.4, 0.1, 0, 0, 2.2, 0, 1.7, 1.0, 0.4, 0.7, 0.2, 0.6, 0, 1.1, 1.1, 1.3, 1.3]
    },
    "LArC": {
        "mode": [102.0, 2.6, 0.7, 1.1, 1.9, 0.9, 0.9, 1.5, 0.4, 2.6, 1.7, 1.7, 1.3, 3.5, 5.9, 2.2, 3.7, 1.5, 1.2, 1.0, 0.5, 0.2, 1.9, 4.8, 2.7, 0.7],
        "low_68": [3.2, 0.5, 0.5, 0, 0.7, 0.4, 0.6, 0, 0.4, 2.3, 0.4, 0.8, 0, 0, 2.1, 1.1, 0.7, 0.4, 0, 0.4, 0.4, 0, 0.7, 1.0, 0.7, 0.5],
        "upp_68": [3.7, 0.8, 0.6, 0, 0.7, 0.7, 0.6, 0, 0.6, 2.1, 0.5, 0.7, 0, 0, 2.5, 1.4, 1.2, 0.5, 0, 0.3, 0.5, 0, 0.8, 1.2, 1.2, 0.7]
    }
}


with open("phII/BEGe_II.json", "w") as json_file:
    json.dump(GERDA_BEGe_phII, json_file, indent=4)

with open("phII/COAX_phII.json", "w") as json_file:
    json.dump(GERDA_COAX_phII, json_file, indent=4)

with open("phIIplus/BEGe_IIplus.json", "w") as json_file:
    json.dump(GERDA_BEGe_phII_plus, json_file, indent=4)

with open("phIIplus/COAX_phIIplus.json", "w") as json_file:
    json.dump(GERDA_COAX_phII_plus, json_file, indent=4)

with open("phIIplus/ICPC_phIIplus.json", "w") as json_file:
    json.dump(GERDA_IC_phII_plus, json_file, indent=4)

# ----------------------------------------------------------------------------------------------
# exposure saved as phII, phII+, phII & phII+
expo_BEGe = [31.46, 21.87, 53.33]
expo_COAX = [28.64, 13.18, 41.82]
expo_ICPC = [0.0, 8.54, 8.54]

# ----------------------------------------------------------------------------------------------
# compute the average among the two periods
def calculate_average_for_det_type(dictII, dictIIplus, det_type):
    # get the correct exposure dict
    if det_type == "BEGe":
        expo_dict = expo_BEGe
    if det_type == "COAX":
        expo_dict = expo_COAX

    # Initialize an empty dictionary to store the averages
    average_dict = {
        "isotope": ["42K", "40K", "60Co", "60Co", "228Ac", "228Ac", "228Ac", "228Ac", "212Bi", "212Pb", "208Tl", "208Tl", "208Tl", "214Pb", "214Pb", "214Pb", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "e+e-", "85Kr", "65Zn"],
        "energy": [1525, 1461, 1333, 1173, 911, 969, 965, 338, 727, 239, 2615, 583, 861, 295, 242, 352, 609, 1765, 1238, 2204, 1378, 2448, 1120, 511, 514, 1125],
        "raw": {
          "mode": [],
          "low_68": [],
          "upp_68": []
        },
        "LArC": {
          "mode": [],
          "low_68": [],
          "upp_68": []
        },
        "LArAC": {
          "mode": [],
          "low_68": [],
          "upp_68": []
        }
    }

    for cut in ["raw", "LArC", "LArAC"]:
        average_dict.update({cut: {}})
        for entry in ["mode", "low_68", "upp_68"]:
            # exposure weightted mean 
            new_entry = [(v1*expo_dict[0] + dictIIplus[cut][entry][idx]*expo_dict[1]) / expo_dict[2] for idx,v1 in enumerate(dictII[cut][entry])]
            # saved in dict
            average_dict[cut].update({entry: new_entry})

    return average_dict

# Calculate the average dictionary
BEGe_average_dictionary = calculate_average_for_det_type(GERDA_BEGe_phII, GERDA_BEGe_phII_plus, "BEGe")
with open("phII_phIIplus/BEGe_II_IIplus.json", "w") as json_file:
    json.dump(BEGe_average_dictionary, json_file, indent=4)

COAX_average_dictionary = calculate_average_for_det_type(GERDA_COAX_phII, GERDA_COAX_phII_plus, "COAX")
with open("phII_phIIplus/COAX_II_IIplus.json", "w") as json_file:
    json.dump(COAX_average_dictionary, json_file, indent=4)

ICPC_average_dictionary = GERDA_IC_phII_plus
with open("phII_phIIplus/ICPC_IIplus.json", "w") as json_file:
    json.dump(ICPC_average_dictionary, json_file, indent=4)

# ----------------------------------------------------------------------------------------------
# calculate total numbers (BEGe + ICPC) for GERDA, through an exposure-weighted mean

# Initialize an empty dictionary to store the averages
BEGe_ICPC_dict = {
    "isotope": ["42K", "40K", "60Co", "60Co", "228Ac", "228Ac", "228Ac", "228Ac", "212Bi", "212Pb", "208Tl", "208Tl", "208Tl", "214Pb", "214Pb", "214Pb", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "214Bi", "e+e-", "85Kr", "65Zn"],
    "energy": [1525, 1461, 1333, 1173, 911, 969, 965, 338, 727, 239, 2615, 583, 861, 295, 242, 352, 609, 1765, 1238, 2204, 1378, 2448, 1120, 511, 514, 1125],
    "raw": {
        "mode": [],
        "low_68": [],
        "upp_68": []
    },
    "LArC": {
        "mode": [],
        "low_68": [],
        "upp_68": []
    },
    "LArAC": {
        "mode": [],
        "low_68": [],
        "upp_68": []
    }
}

for cut in ["raw", "LArC", "LArAC"]:
    BEGe_ICPC_dict.update({cut: {}})
    for entry in ["mode", "low_68", "upp_68"]:
        # exposure weightted mean 
        new_entry = [(v1*expo_BEGe[0] + GERDA_BEGe_phII_plus[cut][entry][idx]*expo_BEGe[1] + GERDA_IC_phII_plus[cut][entry][idx]*expo_ICPC[0] + GERDA_IC_phII_plus[cut][entry][idx]*expo_ICPC[1]) / (expo_BEGe[2] + expo_ICPC[2]) for idx,v1 in enumerate(GERDA_BEGe_phII[cut][entry])]
        # saved in dict
        BEGe_ICPC_dict[cut].update({entry: new_entry})

with open("phII_phIIplus/BEGe_ICPC.json", "w") as json_file:
    json.dump(BEGe_ICPC_dict, json_file, indent=4)