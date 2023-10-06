from legend_plot_style import LEGENDPlotStyle as lps

import os
import json
import re
import numpy as np
import pandas as pd

from legendmeta import JsonDB
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
lps.use('legend_talks')


# function to compute the error of A/B, given A and B errors
def ratio_error(A, A_error, B, B_error):
    return np.sqrt(np.square(A_error/B)+np.square(A*B_error/np.square(B)))


log = False

gamma_energy = {    "K42_1525"  :    1524.7,
                    "K40_1461"  :    1460.8, 
                    "Co60_1332"  :    1332.5, 
                    "Co60_1173"  :    1173.2, 
                    "Ac228_911"  :    911.2, 
                    "Bi212_727"  :    727.3, 
                    "Tl208_2614"  :    2614.5, 
                    "Tl208_583"  :    583.2, 
                    "Tl208_861"  :    860.6, 
                    "Pa234m_1001"  :    1001.0, 
                    "Pb214_352"  :    351.9, 
                    "Pb214_295"  :    295.2, 
                    "Bi214_609"  :    609.3, 
                    "Bi214_1378"  :    1377.7, 
                    "Bi214_1730"  :    1730, 
                    "Bi214_1764"  :    1764.5, 
                    "Bi214_1238"  :    1238.1, 
                    "Bi214_2204"  :    2204.1, 
                    "Bi214_2448"  :    2447.9,
                    "Ac228_1588"  :    1588 }


gamma_isotope = {   "K42_1525"  :    "$^{42}$K",
                    "K40_1461"  :    "$^{40}$K",
                    "Co60_1332"  :    "$^{60}$Co",
                    "Co60_1173"  :    "$^{60}$Co",
                    "Ac228_911"  :    "$^{228}$Ac",
                    "Bi212_727"  :    "$^{212}Bi",
                    "Tl208_2614"  :    "$^{208}$Tl", 
                    "Tl208_583"  :    "$^{208}$Tl", 
                    "Tl208_861"  :    "$^{208}$Tl", 
                    "Pa234m_1001"  :   "$^{234m}$Pa", 
                    "Pb214_352"  :   "$^{214}$Pb", 
                    "Pb214_295"  :    "$^{214}$Pb", 
                    "Bi214_609"  :    "$^{214}$Bi", 
                    "Bi214_1378"  :    "$^{214}$Bi",  
                    "Bi214_1730"  :    "$^{214}$Bi",  
                    "Bi214_1764"  :    "$^{214}$Bi",  
                    "Bi214_1238"  :    "$^{214}$Bi",  
                    "Bi214_2204"  :    "$^{214}$Bi",  
                    "Bi214_2448"  :    "$^{214}$Bi",
                    "Ac228_1588"  :    "$^{228}$Ac" }

periods = "p3_p4"

RAW_path = "/lfs/l1/legend/users/calgaro/taup/legend-gamma-lines-analysis/output/2023-08-21/0_25keV/raw"
AC_path  = "/lfs/l1/legend/users/calgaro/taup/legend-gamma-lines-analysis/output/2023-08-21/0_25keV/LArAC"
C_path   = "/lfs/l1/legend/users/calgaro/taup/legend-gamma-lines-analysis/output/2023-08-21/0_25keV/LArC"

#single_path = f"/lfs/l1/legend/users/calgaro/taup/legend-gamma-lines-analysis/output/2023-08-11/{periods}/potassium_lines/" 
#exposure_path = f"/lfs/l1/legend/users/calgaro/taup/legend-gamma-lines-analysis/livetime_and_exposure/exposure_in_kg_yr_on_no_psd_{periods}.json"

gerda_results = "/lfs/l1/legend/users/morella/gamma-analysis/legend-gamma-lines-analysis/gerda_results/phII_phIIplus/"

# p3 + p4 (messi a mano!!!!!! chaaaaaaaangeeeeee!!!!!)
ICPC_exposure = 8 # kg yr
BEGe_exposure = 2.1 # kg yr

dict_exposure = {"ICPC": ICPC_exposure, "BEGe": BEGe_exposure}

# get channel maps
map_file = os.path.join("/lfs/l1/legend/data/public/prodenv/prod-blind/ref/", "v01.06", "inputs/hardware/configuration/channelmaps")
full_det_map = JsonDB(map_file).on(timestamp="20230601T000000Z")
full_det_map = {k:v for k,v in full_det_map.items() if "location" in v.keys() and "S" not in k} # remove SiPMs and AUX/DUMMY channels

# sort the dictionary keys based on location and position
det_list = sorted(full_det_map.keys(), key=lambda key: (full_det_map[key]['location']['string'], full_det_map[key]['location']['position']))


gamma     = []
cut       = []
det_type  = []
intensity = []
range_min = []
range_max = []

avoid = ["e+e-_Kr85_514", "Pb212_239_Pb214_242", "Ac228_338_Pb214_352", "Ac228_965_Ac228_969", "Bi214_1120_Zn65_1125"]

det_list = ["ICPC"]

for cut_path in [RAW_path, AC_path, C_path]:#, AC_path, C_path]:
    for detector_type in det_list:
            legend = json.load(open("histo.gamma.json"))
            if legend["name_fit"] not in avoid:
                gamma.append("name_fit")
                det_type.append(detector_type)
                if "raw"  in cut_path: cut.append("raw")
                if "LArAC" in cut_path: cut.append("LArAC")
                if "LArC" in cut_path: cut.append("LArC")
                intensity_keys = [k for keys in legend["fit_parameters"]["line"].keys() if "intensity" in keys]:
                for int_key in intensity_keys: 
                    try:
                        intensity.append(legend["fit_parameters"]["line"][int_key]["mode"]/dict_exposure[detector_type])
                        range_min.append(legend["fit_parameters"]["line"][int_key]["low"]/dict_exposure[detector_type])
                        range_max.append(legend["fit_parameters"]["line"][int_key]["high"]/dict_exposure[detector_type])
                    except:
                        intensity.append(legend["fit_parameters"]["line"][int_key]["upper_limit"]/dict_exposure[detector_type])
                        range_min.append(0.0)
                        range_max.append(0.0)
  
df = pd.DataFrame()

df["gamma"]     = gamma
df["det_type"]  = det_type
df["cut"]       = cut
df["intensity"] = intensity
df["range_min"] = range_min
df["range_max"] = range_max

df_2 = df.copy()
df_2["yerr_low"] = df_2.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
df_2["yerr_upp"] = df_2.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
df_2.loc[df_2["intensity"] == df_2["yerr_low"], "yerr_low"] = 0
det_type_order = ["BEGe", "COAX", "ICPC", "PPC"]
df_2['det_type'] = pd.Categorical(df_2['det_type'], categories=det_type_order, ordered=True)
df_2 = df_2.sort_values("det_type")
df = df_2.copy()



legend_color_cut = {"raw" : lps.colors["legend_blue"], "LArAC": lps.colors["legend_darkblue"], "LArC" : lps.colors["legend_orange"]}
gerda_color_cut = {"raw" : lps.colors["legend_orange"], "LArAC": "red"}

label_cut = {"raw": "before LAr", "LArAC": "after LAr", "LArC" : "with LAr"}
legend_dict = {"raw" : "raw", "LArAC": "LAr AC", "LArC" : "LAr C"}

cut_order = ["LArAC", "LArC"]
det_list = ["ICPC"]


# Group the data by gamma and cut
grouped = df.groupby(['gamma', 'cut'])

list_custom_lines = ["K40_1461", "K42_1525", "Bi214_2204", "Tl208_2614"]
list_custom_lines = df.gamma.unique()

marker_size = 5

import re
for g in list_custom_lines:

    # to load gerda data
    isotope_num = g.split("_")[0]
    num = re.findall(r'\d+\.\d+|\d+', isotope_num)[0]
    isotope = isotope_num.split(num)[0]
    first_energy = g.split("_")[1]
    element_to_find = int(first_energy)  
    
    fig, ax = lps.subplots(figsize=(4, 3))
    ax.hlines(1, -0.7, 1.7, color = "black")

    for cut in cut_order: 
        # Plot error bars
        for det_type in det_list:
            condition = (df.gamma == g) & (df.det_type == det_type)

            A_limit_flag = False
            B_limit_flag = False
            lolims       = False

            x = df[condition].det_type.iloc[0]
            
            A = df[condition & (df.cut == "raw")].intensity.iloc[0]
            A_err_low = df[condition & (df.cut == "raw")].yerr_low.iloc[0]
            A_err_upp = df[condition & (df.cut == "raw")].yerr_upp.iloc[0]

            if((A_err_low == 0) & (A_err_upp == 0)): A_limit_flag = True

            B = df[condition & (df.cut == cut)].intensity.iloc[0]
            B_err_low = df[condition & (df.cut == cut)].yerr_low.iloc[0]
            B_err_upp = df[condition & (df.cut == cut)].yerr_upp.iloc[0]

            if((B_err_low == 0) & (B_err_upp == 0)): B_limit_flag = True

            if cut == "LArAC":
                y = B / A
            elif cut == "LArC":
                y = 1 - (B / A)
            yerr_low = ratio_error(B, B_err_low, A, A_err_low)
            yerr_upp = ratio_error(B, B_err_upp, A, A_err_upp)

            if(A_limit_flag | B_limit_flag):
                #lolims = True
                if log == True:
                    yerr_low = 0
                    yerr_upp = 0
                else:
                    yerr_low = 0
                    yerr_upp = 0

            # print(g)
            # print(x)
            # print(y)
            # print(yerr_low)
            # print(yerr_upp)
            # print("\n")

            ax.errorbar(x,
                y,
                yerr=[[yerr_low], [yerr_upp]],
                marker="s", linestyle="", color=legend_color_cut[cut], markersize=marker_size, uplims = lolims
            )
        # Plot error bars

        for det_type in det_list:
            condition = (df.gamma == g) & (df.det_type == det_type)    
            if det_type != "PPC":
                gerda_II = {}
                if det_type == "BEGe":
                    gerda_II = json.load(open(gerda_results + "BEGe_II_IIplus.json"))
                if det_type == "ICPC":
                    gerda_II = json.load(open(gerda_results + "ICPC_IIplus.json"))
                    
                if element_to_find in list(gerda_II["energy"]):
                    index_II = gerda_II["energy"].index(element_to_find)

                    # if gerda_II[cut]["low_68"][index_II] == 0 and gerda_II[cut]["upp_68"][index_II] == 0:

                    #     x = f"{det_type} II"
                    #     y = gerda_II["raw"]["mode"][index_II]/gerda_II["LArAC"]["mode"][index_II] 

                    #     ax.errorbar(x, 
                    #     y, 
                    #     marker="_", 
                    #     uplims = True, 
                    #     yerr = gerda_II[cut]["mode"][index_II], 
                    #     color=legend_color_cut[cut], 
                    #     markersize=marker_size)
                    # else:

                    A_limit_flag = False
                    B_limit_flag = False
                    lolims       = False

                    x = f"{det_type} II"

                    A = gerda_II["raw"]["mode"][index_II]
                    A_err_low = gerda_II["raw"]["low_68"][index_II]
                    A_err_upp = gerda_II["raw"]["upp_68"][index_II]

                    if((A_err_low == 0) & (A_err_upp == 0)): A_limit_flag = True

                    B = gerda_II[cut]["mode"][index_II]
                    B_err_low = gerda_II[cut]["low_68"][index_II]
                    B_err_upp = gerda_II[cut]["upp_68"][index_II]

                    if((B_err_low == 0) & (B_err_upp == 0)): B_limit_flag = True
                    elif(B_err_low == B): B_limit_flag = True

                    if cut == "LArAC":
                        z = B / A
                    elif cut == "LArC":
                        z = 1 - (B / A)

                    zerr_low = ratio_error(B, B_err_low, A, A_err_low)
                    zerr_upp = ratio_error(B, B_err_upp, A, A_err_upp)

                    if(A_limit_flag | B_limit_flag):
                        #lolims = True
                        if log == True:
                            zerr_low = 0
                            zerr_upp = 0
                        else:
                            zerr_low = 0
                            zerr_upp = 0

                    if  z < 0: 
                        print("skipping because of survival < 0") 
                        print(g)
                        print(cut)
                        print(z)
                        continue
                    ax.errorbar(
                        x,
                        z,
                        yerr=[[zerr_low], [zerr_upp]],
                        marker="o", 
                        linestyle="", 
                        color=legend_color_cut[cut], 
                        markersize=marker_size,
                        uplims = lolims
                    )
                    
    ax.grid(axis = "y", color = lps.colors["legend_grey"], linestyle = "--")
    ax.set_xlim((-0.7, 1.7))
    #ax.set_ylim((0, round(max(y,z)*2)))
    if log == True:
        ax.set_ylim((0.01, 2))
        ax.set_yscale("log")
    else: 
        ax.set_ylim((0, 1.2))
    ax.set_xticklabels(["L200", "GERDA"])
    ax.set_ylabel('survival fraction')
    ax.set_legend_logo(position='upper right', logo_type = 'preliminary')
    ax.set_legend_annotation(month="Aug 2023")
    fig.suptitle(f"{gamma_isotope[g]} - {gamma_energy[g]} keV")

    gerda_before  = Line2D([0], [0], marker='o', color = legend_color_cut["raw"], markersize=marker_size)
    gerda_after   = Line2D([0], [0], marker='o', color = legend_color_cut["LArAC"], markersize=marker_size)
    gerda_with    = Line2D([0], [0], marker='o', color = legend_color_cut["LArC"], markersize=marker_size)
    legend_before = Line2D([0], [0], marker='s', color = legend_color_cut["raw"], markersize=marker_size)
    legend_after  = Line2D([0], [0], marker='s', color = legend_color_cut["LArAC"], markersize=marker_size)
    legend_with   = Line2D([0], [0], marker='s', color = legend_color_cut["LArC"], markersize=marker_size)
    
    # if g in ["K40_1461", "Bi214_2204"]:
    #     location = "upper left"
    # else:
    #     location = "center left"
    location = "center left"
    ax.legend([(legend_after, gerda_after), (legend_with, gerda_with)], ["LAr AC", "LAr C"], 
            numpoints=1, 
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = location)
    # if len(cut_order) == 2:
    #     ax.legend([(legend_before, gerda_before), (legend_after, gerda_after)], ["before LAr", "after LAr"], 
    #             numpoints=1, 
    #             handler_map={tuple: HandlerTuple(ndivide=None)},
    #             loc = location)
    # elif len(cut_order) == 3:
    #     ax.legend([(legend_before, gerda_before), (legend_after, gerda_after), (legend_with, gerda_with)], ["before LAr", "after LAr", "with LAr"], 
    #             numpoints=1, 
    #             handler_map={tuple: HandlerTuple(ndivide=None)},
    #             loc = location)
    plt.tight_layout() 
    cwd = os.listdir(os.getcwd())
    if "gamma_survival" not in cwd:
        os.makedirs("gamma_survival")
    #plt.savefig(f"gamma_survival/{g}_survival.pdf")
    plt.savefig(f"gamma_survival/{g}_survival.png")