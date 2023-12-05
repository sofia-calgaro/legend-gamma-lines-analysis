from legend_plot_style import LEGENDPlotStyle as lps
#lps.use('legend_talks')
import os
import uproot
import json
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import sys
sys.path.insert(0, '../src/')
import utils

gamma_isotope = {   "K42_1525"    :    "$^{42}$K",
                    "K40_1461"    :    "$^{40}$K",
                    "Co60_1332"   :    "$^{60}$Co",
                    "Co60_1173"   :    "$^{60}$Co",
                    "Ac228_911"   :    "$^{228}$Ac",
                    "Bi212_727"   :    "$^{212}$Bi",
                    "Tl208_2614"  :    "$^{208}$Tl", 
                    "Tl208_583"   :    "$^{208}$Tl", 
                    "Tl208_861"   :    "$^{208}$Tl", 
                    "Pa234m_1001" :   "$^{234m}$Pa", 
                    "Pb214_352"   :   "$^{214}$Pb", 
                    "Pb214_295"   :    "$^{214}$Pb", 
                    "Bi214_609"   :    "$^{214}$Bi", 
                    "Bi214_1378"  :    "$^{214}$Bi",  
                    "Bi214_1730"  :    "$^{214}$Bi",  
                    "Bi214_1764"  :    "$^{214}$Bi",  
                    "Bi214_1238"  :    "$^{214}$Bi",  
                    "Bi214_2204"  :    "$^{214}$Bi",  
                    "Bi214_2448"  :    "$^{214}$Bi",
                    "Ac228_1588"  :    "$^{228}$Ac",
                    "e+e-"        :    "e$^+$e$^-$",
                    "e+e-_511"    :    "e$^+$e$^-$",
                    "Kr85_514"    :    "$^{85}$Kr",
                    "Pb212_239"   :    "$^{212}$Pb",
                    "Pb214_242"   :    "$^{214}$Pb",
                    "Ac228_338"   :    "$^{228}$Ac",
                    "Pb214_352"   :    "$^{214}$Pb",
                    "Ac228_965"   :    "$^{228}$Ac",
                    "Ac228_969"   :    "$^{228}$Ac",
                    "Bi214_1120"  :    "$^{214}$Bi",
                    "Zn65_1125"   :    "$^{65}$Zn"}

gamma_energy = {    "K42_1525"    :    1524.7,
                    "K40_1461"    :    1460.8, 
                    "Co60_1332"   :    1332.5, 
                    "Co60_1173"   :    1173.2, 
                    "Ac228_911"   :    911.2, 
                    "Bi212_727"   :    727.3, 
                    "Tl208_2614"  :    2614.5, 
                    "Tl208_583"   :    583.2, 
                    "Tl208_861"   :    860.6, 
                    "Pa234m_1001" :    1001.0, 
                    "Pb214_352"   :    351.9, 
                    "Pb214_295"   :    295.2, 
                    "Bi214_609"   :    609.3, 
                    "Bi214_1378"  :    1377.7, 
                    "Bi214_1730"  :    1730, 
                    "Bi214_1764"  :    1764.5, 
                    "Bi214_1238"  :    1238.1, 
                    "Bi214_2204"  :    2204.1, 
                    "Bi214_2448"  :    2447.9,
                    "Ac228_1588"  :    1588,
                    "e+e-"        :    511,
                    "e+e-_511"    :    511.0,
                    "Kr85_514"    :    514.0,
                    "Pb212_239"   :    238.6,
                    "Pb214_242"   :    242.0,
                    "Ac228_338"   :    338.3,
                    "Pb214_352"   :    351.9,
                    "Ac228_965"   :    964.8,
                    "Ac228_969"   :    969.0,
                    "Bi214_1120"  :    1120.3,
                    "Zn65_1125"   :    1125.0}


RAW_path = "/lfs/l1/legend/users/morella/gamma-analysis/output/p03-p04-testing/raw/"
AC_path = "/lfs/l1/legend/users/morella/gamma-analysis/output/p03-p04-testing/LAr_AC/"
C_path = "/lfs/l1/legend/users/morella/gamma-analysis/output/p03-p04-testing/LAr_C/"

RAW_path_canada = "/lfs/l1/legend/users/morella/gamma-analysis/output/vancouver-dataset/raw/"
AC_path_canada = "/lfs/l1/legend/users/morella/gamma-analysis/output/vancouver-dataset/LAr_AC/"
C_path_canada = "/lfs/l1/legend/users/morella/gamma-analysis/output/vancouver-dataset/LAr_C/"

gerda_results = "/lfs/l1/legend/users/morella/gamma-analysis/gamma-testing/attic/gerda_results/phII_phIIplus/"

gamma     = []
cut       = []
det_type  = []
intensity = []
range_min = []
range_max = []
energy_line = []

avoid = ["e+e-_Kr85_514", "Pb212_239_Pb214_242", "Ac228_338_Pb214_352", "Ac228_965_Ac228_969", "Bi214_1120_Zn65_1125"]

dict_exposure = {"ICPC": 9.35, "BEGe": 2.1, "COAX": 1.86, "PPC": 1.42}
dict_exposure.update({"All": sum(dict_exposure.values())})
print("TAUP exposure:", dict_exposure) # CHECK VALUES


for cut_path in [RAW_path, AC_path, C_path]:
    for detector_type in ["ICPC", "BEGe", "PPC", "COAX"]:
        json_files = sorted([file for file in os.listdir(os.path.join(cut_path, detector_type)) if ".json" in file])

        for my_json in json_files: 
            legend = json.load(open(os.path.join(cut_path, detector_type, my_json)))
            gamma_line_name = list(legend.keys())[0]
            det_type.append(detector_type)
            cut.append(cut_path.split("/")[-2])

            if gamma_line_name in avoid: 
                det_type.append(detector_type)
                cut.append(cut_path.split("/")[-2])
                for idx in range(0,2):
                    splits = gamma_line_name.split("_")
                    if len(splits) == 3:
                        if idx==0: name = splits[0]
                        else: name = f"{splits[1]}_{splits[2]}"
                    if len(splits) == 4:
                        if idx==0: name = f"{splits[0]}_{splits[1]}"
                        else: name = f"{splits[2]}_{splits[3]}"
                    if "e+e-" in name: name += "_511"
                    gamma.append(name)
                    energy_line.append(gamma_energy[name])
                    best_fit, cts_L68, cts_U68 = utils.get_results(gamma_line_name, os.path.join(cut_path, detector_type), 1, idx)
                    intensity.append(best_fit/dict_exposure[detector_type])
                    range_min.append(cts_L68/dict_exposure[detector_type])
                    range_max.append(cts_U68/dict_exposure[detector_type])
                continue
            gamma.append(gamma_line_name)
            energy_line.append(gamma_energy[gamma_line_name])
            best_fit, cts_L68, cts_U68 = utils.get_results(gamma_line_name, os.path.join(cut_path, detector_type), 1, 0)
            intensity.append(best_fit/dict_exposure[detector_type])
            range_min.append(cts_L68/dict_exposure[detector_type])
            range_max.append(cts_U68/dict_exposure[detector_type])

df = pd.DataFrame()
df["gamma"]     = gamma
df["det_type"]  = det_type
df["cut"]       = cut
df["energy_line"] = energy_line
df["intensity"] = intensity
df["range_min"] = range_min
df["range_max"] = range_max

df_2 = df.copy()
df_2["yerr_low"] = df_2.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
df_2["yerr_upp"] = df_2.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
df_2.loc[df_2["intensity"] == df_2["yerr_low"], "yerr_low"] = 0
df_2 = df_2.sort_values("det_type")
df = df_2.copy()

# ===================================================================== CANADA
gamma     = []
cut       = []
det_type  = []
energy_line = []
intensity_canada = []
range_min_canada = []
range_max_canada = []
dict_exposure = {}

dict_exposure = {"ICPC": 27.2, "BEGe": 6.3, "COAX": 5.34, "PPC": 4.9}
dict_exposure.update({"All": sum(dict_exposure.values())})
print("vancouver exposure:", dict_exposure)

for cut_path in [RAW_path_canada, AC_path_canada, C_path_canada]:
    for detector_type in ["ICPC", "BEGe", "PPC", "COAX"]:
        if detector_type == "ICPC":
           if "raw" in cut_path:
              cut_path = "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/output_gamma_rates/vancouver_silver_raw"
           if "LAr_C" in cut_path:
              cut_path = "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/output_gamma_rates/vancouver_silver_LArC"
           if "LAr_AC" in cut_path:
              cut_path = "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/output_gamma_rates/vancouver_silver_LArAC"
        else:
           if "raw" in cut_path:
              cut_path = RAW_path_canada
           if "LArC" in cut_path:
              cut_path = C_path_canada
           if "LArAC" in cut_path:
              cut_path = AC_path_canada

        json_files = sorted([file for file in os.listdir(os.path.join(cut_path, detector_type)) if ".json" in file])
        if detector_type=="ICPC": print(os.path.join(cut_path, detector_type), '\n', json_files)

        for my_json in json_files:
            legend = json.load(open(os.path.join(cut_path, detector_type, my_json)))
            gamma_line_name = list(legend.keys())[0]
            det_type.append(detector_type)
            my_cut = ""
            if "raw" in cut_path: my_cut = "raw"
            if "LAr_C" in cut_path or "LArC" in cut_path: my_cut = "LAr_C"
            if "LAr_AC" in cut_path or "LArAC" in cut_path: my_cut = "LAr_AC"
            cut.append(my_cut)

            if gamma_line_name in avoid:
                det_type.append(detector_type)
                cut.append(my_cut)
                for idx in range(0,2):
                    splits = gamma_line_name.split("_")
                    if len(splits) == 3:
                        if idx==0: name = splits[0]
                        else: name = f"{splits[1]}_{splits[2]}"
                    if len(splits) == 4:
                        if idx==0: name = f"{splits[0]}_{splits[1]}"
                        else: name = f"{splits[2]}_{splits[3]}"
                    if "e+e-" in name: name += "_511"
                    gamma.append(name)
                    energy_line.append(gamma_energy[name])
                    best_fit, cts_L68, cts_U68 = utils.get_results(gamma_line_name, os.path.join(cut_path, detector_type), 1, idx)
                    intensity_canada.append(best_fit/dict_exposure[detector_type])
                    range_min_canada.append(cts_L68/dict_exposure[detector_type])
                    range_max_canada.append(cts_U68/dict_exposure[detector_type])
                continue
            gamma.append(gamma_line_name)
            energy_line.append(gamma_energy[gamma_line_name])
            best_fit, cts_L68, cts_U68 = utils.get_results(gamma_line_name, os.path.join(cut_path, detector_type), 1, 0)
            intensity_canada.append(best_fit/dict_exposure[detector_type])
            range_min_canada.append(cts_L68/dict_exposure[detector_type])
            range_max_canada.append(cts_U68/dict_exposure[detector_type])

df_canada = pd.DataFrame()
df_canada["gamma"]     = gamma
df_canada["det_type"]  = det_type
df_canada["cut"]       = cut
df_canada["energy_line"] = energy_line
df_canada["intensity"] = intensity_canada
df_canada["range_min"] = range_min_canada
df_canada["range_max"] = range_max_canada

df_2 = df_canada.copy()
df_2["yerr_low"] = df_2.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
df_2["yerr_upp"] = df_2.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
df_2.loc[df_2["intensity"] == df_2["yerr_low"], "yerr_low"] = 0
df_2 = df_2.sort_values("det_type")
df_canada = df_2.copy()

# =========================================================================== PLOTTING
legend_color_cut = {"raw" : "blue", "LAr_AC": "navy", "LAr_C" : "orange"}
gerda_color_cut = {"raw" : "orange", "LAr_AC": "red"}
label_cut = {"raw": "before LAr", "LArAC": "after LAr", "LArC" : "with LAr"}
legend_dict = {"raw" : "raw", "LAr_AC": "LAr AC", "LAr_C" : "LAr C"}
cut_order = ["raw", "LAr_AC", "LAr_C"]
det_list = ["ICPC", "BEGe", "PPC", "COAX"]
# Group the data by gamma and cut
grouped = df.groupby(['gamma', 'cut'])
list_custom_lines = df.gamma.unique()
list_custom_lines = sorted(list_custom_lines, key=lambda label: gamma_energy[label]) # sorted based on the gamma energy

output_folder = "gamma_rates"

marker_size = 4
import re

df["datasets"] = "TAUP"
df_canada["datasets"] = "Vancouver"
# ========================================================================== plot per gamma, per detector type
"""
for g in list_custom_lines:
    for det_type in det_list:
        fig, ax = lps.subplots(figsize=(5, 4))
        for idx, cut in enumerate(cut_order):

            # TAUP
            condition = (df.gamma == g) & (df.cut == cut) & (df.det_type == det_type)
            yerr_low_zero = (df[condition]["yerr_low"] == 0).all()
            yerr_upp_zero = (df[condition]["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                ax.plot(df[condition].datasets,
                        df[condition].intensity,
                        marker="x", color=legend_color_cut[cut], markersize=marker_size*1.5)
            else:
                ax.errorbar(
                    df[condition].datasets,
                    df[condition].intensity,
                    yerr=[df[condition].yerr_low, df[condition].yerr_upp],
                    marker="s", linestyle="", color=legend_color_cut[cut], markersize=marker_size
                )

            # CANADA
            condition = (df_canada.gamma == g) & (df_canada.cut == cut) & (df_canada.det_type == det_type)
            yerr_low_zero = (df_canada[condition]["yerr_low"] == 0).all()
            yerr_upp_zero = (df_canada[condition]["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                ax.plot(df_canada[condition].datasets,
                         df_canada[condition].intensity,
                         marker="x", color=legend_color_cut[cut], markersize=marker_size*1.5)
            else:
                ax.errorbar(
                    df_canada[condition].datasets,
                    df_canada[condition].intensity,
                    yerr=[df_canada[condition].yerr_low, df_canada[condition].yerr_upp],
                    marker="s", linestyle="", color=legend_color_cut[cut], markersize=marker_size
                )

        ax.set_xlim((-0.7, 1.7))
        ax.set_xticklabels(["TAUP", "Vancouver"])
        ax.set_ylabel('gamma line rate [cts/(kg yr)]')
        #ax.set_legend_logo(position='upper right', logo_type = 'preliminary')
        #ax.set_legend_annotation(month="November 2023")
        ax.grid(axis = "y", color = "grey", linestyle = "--")
        fig.suptitle(f"{gamma_isotope[g]} ({gamma_energy[g]} keV) - {det_type}")

        from matplotlib.lines import Line2D
        from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

        legend_before = Line2D([0], [0], marker='s', color = legend_color_cut["raw"], markersize=marker_size, linewidth=0)
        legend_after  = Line2D([0], [0], marker='s', color = legend_color_cut["LAr_AC"], markersize=marker_size, linewidth=0)
        legend_with   = Line2D([0], [0], marker='s', color = legend_color_cut["LAr_C"], markersize=marker_size, linewidth=0)
        global_mode = Line2D([0], [0], marker='s', color='k', markersize=marker_size)
        upp_lim = Line2D([0], [0], marker='x', color='k', markersize=marker_size*1.5, linewidth=0)

        location = "center left"

        if len(cut_order) == 2:
            ax.legend([(legend_before), (legend_after), (global_mode), (upp_lim)], ["before LAr", "after LAr", "G.M.+-68%", "U90% CI limit"],
                    numpoints=1,
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    loc = location)
        elif len(cut_order) == 3:
            ax.legend([(legend_before), (legend_after), (legend_with), (global_mode), (upp_lim)], ["before LAr", "after LAr", "with LAr", "G.M.+-68%", "U90% CI limit"],
                    numpoints=1,
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    loc = location)
        plt.tight_layout()

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        new_output_folder = os.path.join(output_folder, det_type)
        if not os.path.exists(new_output_folder):
            os.makedirs(new_output_folder)
        plt.savefig(f"{new_output_folder}/{g}_rate_v02.pdf")
        plt.savefig(f"{new_output_folder}/{g}_rate_v02.png")
        
# ========================================================================== plot per gamma, per detector type (TAUP VS VANCOUVER)
df["combined_dataset"] = df["det_type"] + "_" + df["cut"]
df_canada["combined_dataset"] = df_canada["det_type"] + "_" + df_canada["cut"]

for g in list_custom_lines:
    for det_type in det_list:
        fig, ax = lps.subplots(figsize=(5, 4))
        for idx, cut in enumerate(cut_order):

            # TAUP
            condition = (df.gamma == g) & (df.cut == cut) & (df.det_type == det_type)
            yerr_low_zero = (df[condition]["yerr_low"] == 0).all()
            yerr_upp_zero = (df[condition]["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                ax.plot(df[condition].combined_dataset,
                        df[condition].intensity,
                        marker="x", color="dodgerblue", markersize=marker_size*1.5)
            else:
                ax.errorbar(
                    df[condition].combined_dataset,
                    df[condition].intensity,
                    yerr=[df[condition].yerr_low, df[condition].yerr_upp],
                    marker="s", linestyle="", color="dodgerblue", markersize=marker_size
                )

            # CANADA
            condition = (df_canada.gamma == g) & (df_canada.cut == cut) & (df_canada.det_type == det_type)
            yerr_low_zero = (df_canada[condition]["yerr_low"] == 0).all()
            yerr_upp_zero = (df_canada[condition]["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                ax.plot(df_canada[condition].combined_dataset,
                         df_canada[condition].intensity,
                         marker="x", color="orange", markersize=marker_size*1.5)
            else:
                ax.errorbar(
                    df_canada[condition].combined_dataset,
                    df_canada[condition].intensity,
                    yerr=[df_canada[condition].yerr_low, df_canada[condition].yerr_upp],
                    marker="s", linestyle="", color="orange", markersize=marker_size
                )

        #ax.set_xlim((-0.7, 1.7))
        ax.set_xticklabels([f"{det_type}\nraw", f"{det_type}\nLAr AC", f"{det_type}\nLAr C"])
        ax.set_ylabel('gamma line rate [cts/(kg yr)]')
        #ax.set_legend_logo(position='upper right', logo_type = 'preliminary')
        #ax.set_legend_annotation(month="November 2023")
        ax.grid(axis = "y", color = "grey", linestyle = "--")
        fig.suptitle(f"{gamma_isotope[g]} ({gamma_energy[g]} keV) - {det_type}")

        from matplotlib.lines import Line2D
        from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

        legend_before = Line2D([0], [0], marker='s', color = "dodgerblue", markersize=marker_size, linewidth=0)
        legend_after  = Line2D([0], [0], marker='s', color = "orange", markersize=marker_size, linewidth=0)
        global_mode = Line2D([0], [0], marker='s', color='k', markersize=marker_size)
        upp_lim = Line2D([0], [0], marker='x', color='k', markersize=marker_size*1.5, linewidth=0)
        ax.legend([(legend_before), (legend_after), (global_mode), (upp_lim)], ["TAUP", "Vancouver", "G.M.+-68%", "U90% CI limit"],
                numpoints=1,
                handler_map={tuple: HandlerTuple(ndivide=None)},
                loc = location)

        location = "lower left"

        plt.tight_layout()

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        new_output_folder = os.path.join(output_folder, det_type)
        if not os.path.exists(new_output_folder):
            os.makedirs(new_output_folder)
        plt.savefig(f"{new_output_folder}/{g}_rate_v02.pdf")
        plt.savefig(f"{new_output_folder}/{g}_rate_v02.png")

"""
# ================================================================================================ all gammas together, per detector type, per cut
taup_color_cut   = {"raw" : "dodgerblue", "LAr_AC": "dodgerblue", "LAr_C" : "dodgerblue"}
canada_color_cut = {"raw" : "orange", "LAr_AC": "orange", "LAr_C" : "orange"}
output_folder = "gamma_rates_comparison_all_gammas"
my_y_zoom = {"raw": [4e-2, 2.2e2], "LAr_AC": [4e-2, 2.2e2], "LAr_C": [4e-2, 2.2e2]}
print(f"you're going to inspect lines {list_custom_lines}")
for det_type in det_list:
    for idx, cut in enumerate(cut_order):
        fig, ax = lps.subplots(figsize=(10, 7))

        for g in list_custom_lines:
            print(det_type, cut, g)
            # TAUP
            condition = (df.gamma == g) & (df.cut == cut) & (df.det_type == det_type)
            yerr_low_zero = (df[condition]["yerr_low"] == 0).all()
            yerr_upp_zero = (df[condition]["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                ax.plot(df[condition].gamma,
                        df[condition].intensity,
                        marker="x", color=taup_color_cut[cut], markersize=marker_size)
            else:
                ax.errorbar(
                    df[condition].gamma,
                    df[condition].intensity,
                    yerr=[df[condition].yerr_low, df[condition].yerr_upp],
                    marker="s", linestyle="", color=taup_color_cut[cut], markersize=marker_size
                )

            # CANADA
            condition = (df_canada.gamma == g) & (df_canada.cut == cut) & (df_canada.det_type == det_type)
            yerr_low_zero = (df_canada[condition]["yerr_low"] == 0).all()
            yerr_upp_zero = (df_canada[condition]["yerr_upp"] == 0).all()
            
            if yerr_low_zero and yerr_upp_zero:
                ax.plot(df_canada[condition].gamma,
                        df_canada[condition].intensity,
                        marker="x", color=canada_color_cut[cut], markersize=marker_size)
            else:
                ax.errorbar(
                    df_canada[condition].gamma,
                    df_canada[condition].intensity,
                    yerr=[df_canada[condition].yerr_low, df_canada[condition].yerr_upp],
                    marker="s", linestyle="", color=canada_color_cut[cut], markersize=marker_size
                )
            
        #ax.set_xlim((-0.7, 1.7))
        plt.xticks(rotation=90,fontsize=12)
        ax.set_xticklabels(list_custom_lines)
        ax.set_ylabel('gamma line rate [cts/(kg yr)]')
        ax.grid(axis = "y", color = "grey", linestyle = "--")
        fig.suptitle(f"{det_type} - {cut} cut")

        location = "upper right"

        taup_data = Line2D([0], [0], marker='s', color = taup_color_cut[cut], markersize=marker_size, linewidth=0)
        canada_data  = Line2D([0], [0], marker='s', color = canada_color_cut[cut], markersize=marker_size, linewidth=0)
        global_mode = Line2D([0], [0], marker='s', color='k', markersize=marker_size)
        upp_lim = Line2D([0], [0], marker='x', color='k', markersize=marker_size*1.5, linewidth=0)
        plt.ylim(my_y_zoom[cut])
        plt.yscale("log")

        ax.legend([(taup_data), (canada_data), (global_mode), (upp_lim)], ["TAUP", "VANCOUVER", "G.M.+-68%", "U90% CI limit"],
                numpoints=1,
                handler_map={tuple: HandlerTuple(ndivide=None)},
                loc = location)
        plt.tight_layout()

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        plt.savefig(f"{output_folder}/{det_type}_{cut}_rate_v02.pdf")
        plt.savefig(f"{output_folder}/{det_type}_{cut}_rate_v02.png")
    
