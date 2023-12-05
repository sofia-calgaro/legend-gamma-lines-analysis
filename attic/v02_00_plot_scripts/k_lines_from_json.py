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
                    "e+e-_511"    :    511,
                    "Kr85_514"    :    514,
                    "Pb212_239"   :    238.6,
                    "Pb214_242"   :    242,
                    "Ac228_338"   :    338.3,
                    "Pb214_352"   :    351.9,
                    "Ac228_965"   :    964.8,
                    "Ac228_969"   :    969,
                    "Bi214_1120"  :    1120.3,
                    "Zn65_1125"   :    1125}


# THINGS YOU CAN CHANGE - 
dataset = "silver" # silver is giving problems, but I'll fix it if needed - not the priority now
det_inspected = "All" # or "" for BEGe/ICPC/... 

if not os.path.exists("potassium_studies"):
    os.makedirs("potassium_studies")
output_K40 = "potassium_studies/K40"
output_K42 = "potassium_studies/K42"
if not os.path.exists(output_K40):
    os.makedirs(output_K40)
if not os.path.exists(output_K42):
    os.makedirs(output_K42)

RAW_path = f"/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/output_gamma_rates/vancouver_silver_potassium_raw"
AC_path  = f"/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/output_gamma_rates/vancouver_silver_potassium_LArAC"
C_path   = f"/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/output_gamma_rates/vancouver_silver_potassium_LArC"

exposure_path = f"/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/src/settings/exposure_in_kg_yr_on_no_psd.json" # not there!


from legendmeta import JsonDB
map_file = os.path.join("/lfs/l1/legend/data/public/prodenv/prod-blind/ref/", "v02.00", "inputs/hardware/configuration/channelmaps")
full_det_map = JsonDB(map_file).on(timestamp="20230601T000000Z")
full_det_map = {k:v for k,v in full_det_map.items() if "location" in v.keys() and "S" not in k} # remove SiPMs and AUX/DUMMY channels
# sort the dictionary keys based on location and position
det_list = sorted(full_det_map.keys(), key=lambda key: (full_det_map[key]['location']['string'], full_det_map[key]['location']['position']))
exposure_json = json.load(open(exposure_path))
exp_dict = {}
# integrate exposure for each single detector over runs and periods
for det in det_list:
    tmp_exp = 0
    for period in exposure_json.keys():
        for run in exposure_json[period].keys():
            if det not in exposure_json[period][run].keys():
               tmp_exp += 0 # might be off or ac
            else:
               tmp_exp += exposure_json[period][run][det]
    exp_dict[det] = tmp_exp
exp_dict = {k:v for k,v in exp_dict.items() if v != 0} # remove dets for which the exposure is null\n",
icpc = 0
coax = 0
bege = 0
ppc = 0
for k in exp_dict.keys():
    if k[0] == "B":
        bege+=exp_dict[k]
    if k[0] == "C":
        coax+=exp_dict[k]
    if k[0] == "P":
        ppc+=exp_dict[k]
    if k[0] == "V":
        icpc+=exp_dict[k]
print (f"bege:", bege)
print (f"coax:", coax)
print (f"icpc:", icpc)
print(f"ppc:", ppc)
#print(exp_dict)


detector4  = []
intensity4 = []
range_min4 = []
range_max4 = []
gamma4     = []
energy_line4 = []
cut4       = []
for cut_path in [RAW_path, AC_path, C_path]:
    for det in os.listdir(os.path.join(cut_path, "single")):
        json_files = sorted([file for file in os.listdir(os.path.join(cut_path, "single", det)) if ".json" in file])

        for my_json in json_files:
            legend = json.load(open(os.path.join(cut_path, "single", det, my_json)))
            gamma_line_name = list(legend.keys())[0]
            detector4.append(det)
            cut4.append(cut_path.split("_")[-1])
            gamma4.append(gamma_line_name)
            energy_line4.append(gamma_energy[gamma_line_name])
            best_fit, cts_L68, cts_U68 = utils.get_results(gamma_line_name, os.path.join(cut_path, "single", det), 1, 0)
            intensity4.append(best_fit/exp_dict[det])
            range_min4.append(cts_L68/exp_dict[det])
            range_max4.append(cts_U68/exp_dict[det])

print(len(gamma4))
print(len(cut4))
print(len(detector4))
print(len(intensity4))
print(len(range_min4))
print(len(range_max4))
single_df = pd.DataFrame()
single_df["gamma"] = gamma4
single_df["cut"] = cut4
single_df["detector"] = detector4
single_df["intensity"] = intensity4
single_df["energy_line"] = energy_line4
single_df["range_min"] = range_min4
single_df["range_max"] = range_max4

from legendmeta import JsonDB
map_file = os.path.join("/lfs/l1/legend/data/public/prodenv/prod-blind/ref/v02.00", "inputs/dataprod/config")
full_status_map = JsonDB(map_file).on(
    timestamp="20230701T000000Z", system="geds"
)["analysis"]
status_map = {key: value for key, value in full_status_map.items() if value.get('usability') in ["off", "ac", "no_psd"]}
all_detectors = [det for det in status_map.keys() if "location" in status_map[det].keys() and "S" not in det]# remove SiPMs and AUX/DUMMY channels
map_file = os.path.join("/lfs/l1/legend/data/public/prodenv/prod-blind/ref/v02.00", "inputs/hardware/configuration/channelmaps")
channel_map = JsonDB(map_file).on(timestamp="20230701T000000Z")

def get_position(detector_key):
    return channel_map[detector_key]['location']['position']
def get_string(detector_key):
    return channel_map[detector_key]['location']['string']
def get_expo(detector_key):
    return exp_dict[detector_key]
def get_det_type(detector_key):
    if "C0" in detector_key: i_det_type = "COAX"
    if "B0" in detector_key: i_det_type = "BEGe"
    if "P0" in detector_key: i_det_type = "PPC"
    if "V0" in detector_key: i_det_type = "ICPC"
    return i_det_type

single_df['string'] = single_df['detector'].apply(get_string)
single_df['position'] = single_df['detector'].apply(get_position)
single_df['expo'] = single_df['detector'].apply(get_expo)
single_df['name'] = "s"+single_df['string'].astype(str)+"-p"+single_df['position'].astype(str)+"-"+single_df['detector']
single_df['det_type'] = single_df['detector'].apply(get_det_type)

color_cut = {"raw" : "royalblue", "LArAC":  "navy", "LArC" : "orange"}
legend_dict = {"raw" : "raw", "LArAC": "LAr AC", "LArC" : "LAr C"}
cut_order = ["LArC", "LArAC", "raw"]
peaks = ["K40_1461", "K42_1525"]

# ======================================================================================================================= plot rate per each detector 
for peak in peaks:
    fix, axes = plt.subplots(figsize=(18,8))
    axes.set_legend_logo(position='upper right')
    for cut in cut_order:
        df = single_df[(single_df.gamma == peak) & (single_df.cut == cut)].copy()
        df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0
        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df = df.sort_values("detector")
        

        for det in df.detector.values:
            yerr_low_zero = (df[df["detector"] == det]["yerr_low"] == 0).all()
            yerr_upp_zero = (df[df["detector"] == det]["yerr_upp"] == 0).all()
            
            if det == df.detector.values[0]:    
                if yerr_low_zero and yerr_upp_zero:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["intensity"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[cut], markersize=5, label=legend_dict[cut]
                    )
                else:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["intensity"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        marker=".", linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                    )
            else:
                if yerr_low_zero and yerr_upp_zero:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["intensity"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, capsize=3
                        marker="x", color=color_cut[cut], markersize=5
                    )
                else:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["intensity"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        marker=".", linestyle="", color=color_cut[cut], markersize=5, capsize=3
                    )

    # delimit strings
    first_dets_of_a_string = []
    if dataset == "silver":
        first_dets_of_a_string = [0, 9, 9+8, 9+8+7, 9+8+7+7, 3+9+8+6+7+3, 3+9+8+6+7+3+7, 3+9+8+6+7+3+7+8, 5+9+8+6+7+3+7+8+7, 7+9+8+6+7+3+7+8+7+11]
    else:
        first_dets_of_a_string = []
    for string in first_dets_of_a_string:
        if string == 0:  continue
        if string == 52: continue
        plt.axvline(string-0.5, linestyle='--', color='gray')
    if peak == "K40_1461": y_text = 300
    if peak == "K42_1525": y_text = 410
    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    plt.xticks(rotation=90, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("rate [counts/(kg yr)]")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v], 
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'], 
            numpoints=1, 
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    plt.tight_layout()
    if peak == "K40_1461": plt.ylim(-5,y_text+10-100)
    if peak == "K42_1525": plt.ylim(-5,y_text+15-100)
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_single.png")
        plt.savefig(f"{output_K40}/K40_single.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_single.png")
        plt.savefig(f"{output_K42}/K42_single.pdf")


# ======================================================================================================================= plot rate per each detector VS external surface 
my_strings = [1,2,3,4,5,7,8,9,10,11]
from legendmeta import LegendMetadata
from datetime import datetime
lmeta = LegendMetadata()
map_file = os.path.join("/lfs/l1/legend/data/public/prodenv/prod-blind/ref/", "v02.00", "inputs/hardware/configuration/channelmaps")
first_key = '20230311T235840Z' # start of p03
chmap  = JsonDB(map_file).on(timestamp=first_key)
channels = [chn for chn in chmap.keys() if chmap[chn]['system']=='geds'] # select HPGe channels
strings = np.array([chmap[chn]['location']['string'] for chn in channels])
"""
# get external surfaces
ge_surface = {}
for stirng in my_strings:
    for chn,val in lmeta.channelmap(datetime.now()).items():
        if val.system == "geds":
            if chmap[chn]['location']['string'] == string:
                ge_surface.update({chn: {"surface": }})
"""
def get_surface(detector_key):
    my_chn = lmeta.channelmap(datetime.now())[detector_key]
    r_mm = my_chn.geometry.radius_in_mm
    h_mm = my_chn.geometry.height_in_mm
    return 2*np.pi*r_mm*h_mm

single_df['ge_surface_mm_2'] = single_df['detector'].apply(get_surface)

for peak in peaks:
    fix, axes = plt.subplots(figsize=(18,8))
    axes.set_legend_logo(position='upper right')
    for cut in cut_order:
        df = single_df[(single_df.gamma == peak) & (single_df.cut == cut)].copy()
        df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0
        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the external surface
        df = df.sort_values("ge_surface_mm_2")


        for det in df.detector.values:
            yerr_low_zero = (df[df["detector"] == det]["yerr_low"] == 0).all()
            yerr_upp_zero = (df[df["detector"] == det]["yerr_upp"] == 0).all()

            if det == df.detector.values[0]:
                if yerr_low_zero and yerr_upp_zero:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["intensity"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[cut], markersize=5, label=legend_dict[cut]
                    )
                else:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["intensity"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        marker=".", linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                    )
            else:
                if yerr_low_zero and yerr_upp_zero:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["intensity"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, capsize=3
                        marker="x", color=color_cut[cut], markersize=5
                    )
                else:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["intensity"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        marker=".", linestyle="", color=color_cut[cut], markersize=5, capsize=3
                    )

    if peak == "K40_1461": y_text = 300
    if peak == "K42_1525": y_text = 410
    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    plt.xticks(rotation=90, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("rate [counts/(kg yr)]")
    plt.xlabel("detectors sorted by external surface dimension")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v],
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'],
            numpoints=1,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    plt.tight_layout()
    if peak == "K40_1461": plt.ylim(-5,y_text+10-100)
    if peak == "K42_1525": plt.ylim(-5,y_text+15-100)
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_surface.png")
        plt.savefig(f"{output_K40}/K40_surface.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_surface.png")
        plt.savefig(f"{output_K42}/K42_surface.pdf")

# ======================================================================================================================= plot survival fraction per each detector
# survival fractions for each detector
single_df_raw = single_df.copy()[single_df.copy()['cut'] == "raw"]
single_df_LArAC = single_df.copy()[single_df.copy()['cut'] == "LArAC"]
single_df_LArC = single_df.copy()[single_df.copy()['cut'] == "LArC"]

merged_df_LArAC = pd.merge(single_df_raw, single_df_LArAC, on=["detector", "gamma", "string", "position", "expo", "name", "det_type"], suffixes=('_raw', '_LArAC'))
merged_df_LArC = pd.merge(single_df_raw, single_df_LArC, on=["detector", "gamma", "string", "position", "expo", "name", "det_type"], suffixes=('_raw', '_LArC'))

merged_df_LArAC['survival_fr'] = merged_df_LArAC['intensity_LArAC'] / merged_df_LArAC['intensity_raw']
merged_df_LArC['survival_fr'] = merged_df_LArC['intensity_LArC'] / merged_df_LArC['intensity_raw']
merged_df_LArC['survival_fr'] = 1 - merged_df_LArC['survival_fr']

# remove duplicated columns
merged_df_LArAC = merged_df_LArAC.drop(columns=['cut_raw', 'range_min_LArAC', 'range_max_LArAC'])
merged_df_LArC = merged_df_LArC.drop(columns=['cut_raw', 'range_min_LArC', 'range_max_LArC'])

# rename columns
merged_df_LArAC = merged_df_LArAC.rename(columns={'cut_LArAC': "cut"})
merged_df_LArC = merged_df_LArC.rename(columns={'cut_LArC': "cut"})


sf_cut_order = ["LArC", "LArAC"]
for peak in peaks:
    fix, axes = plt.subplots(figsize=(18,8))
    axes.set_legend_logo(position='upper right')
    for cut in sf_cut_order:
        if cut == "LArC":
            df = merged_df_LArC[(merged_df_LArC.gamma == peak)].copy()
        else:
            df = merged_df_LArAC[(merged_df_LArAC.gamma == peak)].copy()
        #df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        #df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        #df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0
        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df = df.sort_values("detector")
        print(df.to_markdown())

        # remove me in the future
        df["yerr_low"] = 0
        df["yerr_upp"] = 0

        for det in df.detector.values:
            if det == df.detector.values[0]:
                if df[df.detector==det]["yerr_low"].values[0]!=0 and df[df.detector==det]["yerr_upp"].values[0]!=0 and df[df.detector==det]["survival_fr"].values[0]!=0:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["survival_fr"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        marker=".", linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                    )
                else:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["survival_fr"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[cut], markersize=5, label=legend_dict[cut]
                    )
            else:
                if df[df.detector==det]["yerr_low"].values[0]!=0 and df[df.detector==det]["yerr_upp"].values[0]!=0 and df[df.detector==det]["survival_fr"].values[0]!=0:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["survival_fr"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        marker=".", linestyle="", color=color_cut[cut], markersize=5, capsize=3
                    )
                else:
                    plt.errorbar(
                        df[df.detector==det]["name"],
                        df[df.detector==det]["survival_fr"],
                        yerr=[df[df.detector==det]["yerr_low"], df[df.detector==det]["yerr_upp"]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, capsize=3
                        marker="x", color=color_cut[cut], markersize=5
                    )

    # delimit strings
    first_dets_of_a_string = []
    if dataset == "silver":
        first_dets_of_a_string = [0, 9, 9+8, 9+8+7, 9+8+7+7, 3+9+8+6+7+3, 3+9+8+6+7+3+7, 3+9+8+6+7+3+7+8, 5+9+8+6+7+3+7+8+7, 7+9+8+6+7+3+7+8+7+11]
    else:
        first_dets_of_a_string = []
    for string in first_dets_of_a_string:
        if string == 0:  continue
        if string == 52: continue
        plt.axvline(string-0.5, linestyle='--', color='gray')
    #if peak == "K40_1461": y_text = 300
    #if peak == "K42_1525": y_text = 410

    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    plt.xticks(rotation=90, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("survival fraction")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v],
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'],
            numpoints=1,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    plt.tight_layout()
    #if peak == "K40_1461": plt.ylim(-5,250)
    #if peak == "K42_1525": plt.ylim(-5,350)
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_single_SF.png")
        plt.savefig(f"{output_K40}/K40_single_SF.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_single_SF.png")
        plt.savefig(f"{output_K42}/K42_single_SF.pdf")


# ======================================================================================================================= plot VS position in the array
with open('LGND_200_Baseline.json', 'r') as json_file:
    LGND_200_BSLN = json.load(json_file)

crys_assem = [k for k in LGND_200_BSLN.keys() if "Crys" in k]
crys_assem_dets = [d.split('CrysAssem')[-1] for d in crys_assem]
rod_heights = {crys_assem_dets[idx]: LGND_200_BSLN[c]['parameters']['rodheight'] for idx,c in enumerate(crys_assem)}

strings = [channel_map[k]["location"]["string"] for k,v in channel_map.items() if "location" in v and "string" in v["location"]]
strings = list(set(strings))
new_rod_heights = {}
for string in strings:
    position = [channel_map[k]["location"]["position"] for k,v in channel_map.items() if "location" in v and "string" in v["location"] and v["location"]["string"]==string]
    dets = [k for k,v in channel_map.items() if "location" in v and "string" in v["location"] and v["location"]["string"]==string]
    
    current_position = 0
    for det in dets:
        current_position += float(rod_heights[det])
        new_rod_heights.update({det: current_position})
print(new_rod_heights)

def get_z_pos(detector_key):
    return new_rod_heights[detector_key]

single_df['z_pos_mm'] = single_df['detector'].apply(get_z_pos)
#single_df[(single_df.string==1) & (single_df.cut=="raw") & (single_df.gamma=="K40_1461")]
all_rods = [v for v in new_rod_heights.values()] # full array
print("min:", min(all_rods))
print("max:", max(all_rods))
plt.hist(all_rods, bins=60)  
plt.axvline(x=130)
plt.axvline(x=300)
plt.axvline(x=450)
plt.axvline(x=730)
plt.xlabel('z position (mm)')
plt.ylabel('frequency')
plt.title('Full array - distribution of z positions')
plt.savefig(f"potassium_studies/z_pos.png")

def determine_z_pos(z_pos_mm):
    if z_pos_mm <= 130:
        return 'top'
    if z_pos_mm > 130 and z_pos_mm <= 300:
        return 'mid-top'
    if z_pos_mm > 300 and z_pos_mm <= 450:
        return 'middle'
    if z_pos_mm > 450 and z_pos_mm <= 730:
        return 'mid-bottom'
    if z_pos_mm >= 730:
        return 'bottom'

single_df['z_pos'] = single_df['z_pos_mm'].apply(determine_z_pos)
for k,z_pos_mm in new_rod_heights.items():
    if z_pos_mm >= 730:
        print(k, channel_map[k]["location"]["string"], channel_map[k]["location"]["position"])
    
####################################################################################################
for peak in peaks:
    fix, axes = plt.subplots(figsize=(10,5))
    axes.set_legend_logo(position='upper right')
    for my_cut in ["raw", "LArC", "LArAC"]:
        df = single_df[(single_df.gamma == peak) & (single_df.cut == my_cut)].copy()

        df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0

        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df = df.sort_values("detector")

        for string_det in list(df.z_pos.unique()):
            df_string = df[df["z_pos"]==string_det]
            # expo-weighted intensity
            intensity=0
            low_err=0
            upp_err=0
            for index, row in df_string.iterrows():
                intensity += row['intensity'] * row['expo']
                low_err += row['yerr_low'] * row['expo']
                upp_err += row['yerr_upp'] * row['expo']
            intensity /= df_string["expo"].sum()
            low_err /= df_string["expo"].sum()
            upp_err /= df_string["expo"].sum()

            my_label = ""
            if string_det == 1: my_label = legend_dict[cut]
            
            yerr_low_zero = (df_string["yerr_low"] == 0).all()
            yerr_upp_zero = (df_string["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                plt.errorbar(
                        str(string_det),
                        intensity,
                        yerr=[[low_err], [upp_err]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[my_cut], markersize=5, label=my_label
                )
            else: 
                plt.errorbar(
                        str(string_det),
                        intensity,
                        yerr=[[low_err], [upp_err]], 
                        marker=".", linestyle="", color=color_cut[my_cut], markersize=5, label=my_label, capsize=3
                )
            
    if peak == "K40_1461": 
        y_text = 170
    if peak == "K42_1525":
        y_text = 360
    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    #plt.xticks(rotation=90, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("exposure weighted rate [counts/(kg yr)]")
    plt.xlabel("position in the array")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v],
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'],
            numpoints=1,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    if peak == "K40_1461": plt.ylim(-5,y_text+10)
    if peak == "K42_1525": plt.ylim(-5,y_text+15)
    #plt.show()
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_position.png")
        plt.savefig(f"{output_K40}/K40_position.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_position.png")
        plt.savefig(f"{output_K42}/K42_position.pdf")



# ======================================================================================================================= plot VS LAr volume in MSs
mss = [
    [956  , 51,  109   ,  188.79, -25],    # 1
    [960  , 48,  198.22,   95.44, -27],    # 2
    [652  , 44,  220   ,   0    , 127],    # 3
    [898  , 48,  198.22,  -95.44,   4],    # 4
    [543  , 54,  109   , -188.79, 181.5],  # 5
    [795  , 51, -109   , -188.79,  55.5],  # 7
    [959.5, 48, -198.22,  -95.44, -26.75], # 8
    [966.5, 44, -220   ,    0   , -30.25], # 9
    [951.5, 48, -198.22,   95.44, -22.75], # 10
    [913  , 51, -109   ,  188.79,  -3.5],  # 11
]
my_strings = [1,2,3,4,5,7,8,9,10,11]
from legendmeta import LegendMetadata
from datetime import datetime
lmeta = LegendMetadata()
map_file = os.path.join("/lfs/l1/legend/data/public/prodenv/prod-blind/ref/", "v02.00", "inputs/hardware/configuration/channelmaps")
first_key = '20230311T235840Z' # start of p03
chmap  = JsonDB(map_file).on(timestamp=first_key)
channels = [chn for chn in chmap.keys() if chmap[chn]['system']=='geds'] # select HPGe channels
strings = np.array([chmap[chn]['location']['string'] for chn in channels])

# get LAr volume of each MS removing from MS volume the Ge volume
ar_volumes = []
ar_volumes_text = []
for idx,ms_info in enumerate(mss):
    h_ms = ms_info[0]
    r_ms = ms_info[1]
    vol_ms = (np.pi*r_ms**2)*h_ms
    vol_string = 0
    string = my_strings[idx]
    for chn,val in lmeta.channelmap(datetime.now()).items():
        if val.system == "geds":
            if chmap[chn]['location']['string'] == string:
                vol_string += (np.pi*val.geometry.radius_in_mm**2)*val.geometry.height_in_mm
    ar_volumes.append(round((vol_ms-vol_string)*1e-6, 2))
    ar_volumes_text.append(f"{round((vol_ms-vol_string)*1e-6, 2)}\nstr.{string}") # m^3

print("Ar volume:", ar_volumes)
print("These are my strings:", my_strings)

yxz = zip(ar_volumes, ar_volumes_text, my_strings)
yxz = sorted(yxz)
sorted_ar_volumes = [y for y,x,z in yxz]
sorted_ar_volumes_text = [x for y,x,z in yxz]
sorted_my_strings = [z for y,x,z in yxz]

print("Ar volume:", sorted_ar_volumes)
print("These are my strings:", sorted_my_strings)

for peak in peaks:
    fix, axes = plt.subplots(figsize=(10,5))
    axes.set_legend_logo(position='upper right')
    for my_cut in ["raw", "LArC", "LArAC"]:
        df = single_df[(single_df.gamma == peak) & (single_df.cut == my_cut)].copy()

        df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0

        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df = df.sort_values("detector")

        for string_det in sorted_my_strings:
            my_str = sorted_my_strings.index(string_det)
            df_string = df[df["string"]==string_det]
            # expo-weighted intensity
            intensity=0
            low_err=0
            upp_err=0
            for index, row in df_string.iterrows():
                intensity += row['intensity'] * row['expo']
                low_err += row['yerr_low'] * row['expo']
                upp_err += row['yerr_upp'] * row['expo']
            intensity /= df_string["expo"].sum()
            low_err /= df_string["expo"].sum()
            upp_err /= df_string["expo"].sum()

            my_label = ""
            if string_det == 1: my_label = legend_dict[cut]

            yerr_low_zero = (df_string["yerr_low"] == 0).all()
            yerr_upp_zero = (df_string["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                plt.errorbar(
                        str(sorted_ar_volumes_text[my_str]),
                        intensity,
                        yerr=[[low_err], [upp_err]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[my_cut], markersize=5, label=my_label
                )

            else:
                plt.errorbar(
                        str(sorted_ar_volumes_text[my_str]),
                        intensity,
                        yerr=[[low_err], [upp_err]],
                        marker=".", linestyle="", color=color_cut[my_cut], markersize=5, label=my_label, capsize=3
                )


    if peak == "K40_1461":
        y_text = 170
    if peak == "K42_1525":
        y_text = 360
    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    plt.xticks(rotation=30, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("exposure weighted rate [counts/(kg yr)]")
    plt.xlabel("LAr volume (m^3)")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v],
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'],
            numpoints=1,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    if peak == "K40_1461": plt.ylim(-5,y_text+10)
    if peak == "K42_1525": plt.ylim(-5,y_text+15)
    plt.tight_layout()
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_LAr_volume.png")
        plt.savefig(f"{output_K40}/K40_LAr_volume.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_LAr_volume.png")
        plt.savefig(f"{output_K42}/K42_LAr_volume.pdf")

# ======================================================================================================================= plot VS string length
my_strings = [1,2,3,4,5,7,8,9,10,11]
crys_assem = [k for k in LGND_200_BSLN.keys() if "Crys" in k]
crys_assem_dets = [d.split('CrysAssem')[-1] for d in crys_assem]
rod_heights = {crys_assem_dets[idx]: LGND_200_BSLN[c]['parameters']['rodheight'] for idx,c in enumerate(crys_assem)}

string_length = []
string_length_text = []
for string in my_strings:
    position = [channel_map[k]["location"]["position"] for k,v in channel_map.items() if "location" in v and "string" in v["location"] and v["location"]["string"]==string]
    dets = [k for k,v in channel_map.items() if "location" in v and "string" in v["location"] and v["location"]["string"]==string]
    
    tot_length = 0
    for det in dets:
        tot_length += float(rod_heights[det])
    string_length.append(tot_length)
    string_length_text.append(f"{tot_length}\nstr.{string}")

print("String lengths:", string_length)
print("These are my strings:", my_strings)

xyz = sorted(zip(string_length, my_strings, string_length_text))
sorted_string_length = [x for x,y,z in xyz]
sorted_my_strings = [y for x,y,z in xyz]
sorted_string_length_text = [z for x,y,z in xyz]

for peak in peaks:
    fix, axes = plt.subplots(figsize=(10,5))
    axes.set_legend_logo(position='upper right')
    for my_cut in ["raw", "LArC", "LArAC"]:
        df = single_df[(single_df.gamma == peak) & (single_df.cut == my_cut)].copy()

        df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0

        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df = df.sort_values("detector")

        for string_det in sorted_my_strings:
            my_str = sorted_my_strings.index(string_det)
            df_string = df[df["string"]==string_det]
            # expo-weighted intensity
            intensity=0
            low_err=0
            upp_err=0
            for index, row in df_string.iterrows():
                intensity += row['intensity'] * row['expo']
                low_err += row['yerr_low'] * row['expo']
                upp_err += row['yerr_upp'] * row['expo']
            intensity /= df_string["expo"].sum()
            low_err /= df_string["expo"].sum()
            upp_err /= df_string["expo"].sum()

            my_label = ""
            if string_det == 1: my_label = legend_dict[cut]

            yerr_low_zero = (df_string["yerr_low"] == 0).all()
            yerr_upp_zero = (df_string["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                plt.errorbar(
                        str(sorted_string_length_text[my_str]),
                        intensity,
                        yerr=[[low_err], [upp_err]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[my_cut], markersize=5, label=my_label
                )

            else:
                plt.errorbar(
                        str(sorted_string_length_text[my_str]),
                        intensity,
                        yerr=[[low_err], [upp_err]],
                        marker=".", linestyle="", color=color_cut[my_cut], markersize=5, label=my_label, capsize=3
                )


    if peak == "K40_1461":
        y_text = 170
    if peak == "K42_1525":
        y_text = 360
    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    plt.xticks(rotation=30, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("exposure weighted rate [counts/(kg yr)]")
    plt.xlabel("Total string length (mm)")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v],
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'],
            numpoints=1,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    if peak == "K40_1461": plt.ylim(-5,y_text+10)
    if peak == "K42_1525": plt.ylim(-5,y_text+15)
    plt.tight_layout()
    #plt.show()
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_string_length.png")
        plt.savefig(f"{output_K40}/K40_string_length.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_string_length.png")
        plt.savefig(f"{output_K42}/K42_string_length.pdf")
exit()    
# ======================================================================================================================= plot VS position, separately for each cut    
"""
for peak in peaks:
    for my_cut in ["raw", "LArC", "LArAC"]:
        fix, axes = plt.subplots(figsize=(10,5))
        axes.set_legend_logo(position='upper right')
        df = single_df[(single_df.gamma == peak) & (single_df.cut == my_cut)].copy()

        df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0

        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df = df.sort_values("detector")

        for string_det in df.z_pos:
            df_string = df[df["z_pos"]==string_det]
            # expo-weighted intensity
            intensity=0
            low_err=0
            upp_err=0
            for index, row in df_string.iterrows():
                intensity += row['intensity'] * row['expo']
                low_err += row['yerr_low'] * row['expo']
                upp_err += row['yerr_upp'] * row['expo']
            intensity /= df_string["expo"].sum()
            low_err /= df_string["expo"].sum()
            upp_err /= df_string["expo"].sum()
            
            if string_det == 1:
                plt.errorbar(
                    str(string_det),
                    intensity,
                    yerr=[[low_err], [upp_err]],
                    marker=".", linestyle="", color=color_cut[my_cut], markersize=5, label=legend_dict[my_cut], capsize=3
                )
            else:
                plt.errorbar(
                    str(string_det),
                    intensity,
                    yerr=[[low_err], [upp_err]],
                    marker=".", linestyle="", color=color_cut[my_cut], markersize=5, capsize=3
                )

        if peak == "K40_1461": 
            if my_cut == "raw": y_text = 150
            if my_cut == "LArC": y_text = 20
            if my_cut == "LArAC": y_text = 140
        if peak == "K42_1525":
            if my_cut == "raw": y_text = 380
            if my_cut == "LArC": y_text = 280
            if my_cut == "LArAC": y_text = 80

        # Customize x-axis tick labels color based on the starting letter
        tick_labels = plt.xticks()[1]
        #axes.set_legend_annotation()
        #plt.xticks(rotation=90, ha='center')
        if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
        if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
        plt.ylabel("exposure weighted rate [counts/(kg yr)]")
        plt.xlabel("position in the array")

        marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
        marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
        marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
        marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
        if my_cut == "raw": 
            axes.legend([(marker_proxy_dot_raw)], 
                    ["raw"], 
                    numpoints=1, 
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    loc = "upper left")
        if my_cut == "LArC":
            axes.legend([(marker_proxy_dot_LArC)], 
                    ["LAr C"], 
                    numpoints=1, 
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    loc = "upper left")
        if my_cut == "LArAC": 
            axes.legend([(marker_proxy_dot_LArAC)], 
                    ["LAr AC"], 
                    numpoints=1, 
                    handler_map={tuple: HandlerTuple(ndivide=None)},
                    loc = "upper left")

        if peak == "K40_1461": plt.ylim(0,y_text+10)
        if peak == "K42_1525": plt.ylim(0,y_text+15)
        #plt.show()
        if peak == "K40_1461":
            plt.savefig(f"{output_K40}/K40_position_{my_cut}.png")
            plt.savefig(f"{output_K40}/K40_position_{my_cut}.pdf")
        else:
            plt.savefig(f"{output_K42}/K42_position_{my_cut}.png")
            plt.savefig(f"{output_K42}/K42_position_{my_cut}.pdf")
"""

# ======================================================================================================================= plot VS string
for peak in peaks:
    fix, axes = plt.subplots(figsize=(10,5))
    axes.set_legend_logo(position='upper right')
    for my_cut in ["raw", "LArC", "LArAC"]:
        df = single_df[(single_df.gamma == peak) & (single_df.cut == my_cut)].copy()

        df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0

        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df = df.sort_values("detector")

        for string_det in strings:
            df_string = df[df["string"]==string_det]
            # expo-weighted intensity
            intensity=0
            low_err=0
            upp_err=0
            for index, row in df_string.iterrows():
                intensity += row['intensity'] * row['expo']
                low_err += row['yerr_low'] * row['expo']
                upp_err += row['yerr_upp'] * row['expo']
            intensity /= df_string["expo"].sum()
            low_err /= df_string["expo"].sum()
            upp_err /= df_string["expo"].sum()

            my_label = ""
            if string_det == 1: my_label = legend_dict[cut]
            
            yerr_low_zero = (df_string["yerr_low"] == 0).all()
            yerr_upp_zero = (df_string["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                plt.errorbar(
                        str(string_det),
                        intensity,
                        yerr=[[low_err], [upp_err]], 
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[my_cut], markersize=5, label=my_label
                )

            else:
                plt.errorbar(
                        str(string_det),
                        intensity,
                        yerr=[[low_err], [upp_err]],
                        marker=".", linestyle="", color=color_cut[my_cut], markersize=5, label=my_label, capsize=3
                )
            
            
    if peak == "K40_1461": 
        y_text = 170
    if peak == "K42_1525":
        y_text = 360
    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    #plt.xticks(rotation=90, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("exposure weighted rate [counts/(kg yr)]")
    plt.xlabel("string")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v],
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'],
            numpoints=1,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    if peak == "K40_1461": plt.ylim(-5,y_text+10)
    if peak == "K42_1525": plt.ylim(-5,y_text+15)
    #plt.show()
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_string.png")
        plt.savefig(f"{output_K40}/K40_string.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_string.png")
        plt.savefig(f"{output_K42}/K42_string.pdf")



# ======================================================================================================================= plot VS detector type
for peak in peaks:
    fix, axes = plt.subplots(figsize=(10,5))
    axes.set_legend_logo(position='upper right')
    for my_cut in ["raw", "LArC", "LArAC"]:
        df = single_df[(single_df.gamma == peak) & (single_df.cut == my_cut)].copy()

        df["yerr_low"] = df.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df["yerr_upp"] = df.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df.loc[df["intensity"] == df["yerr_low"], "yerr_low"] = 0

        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df["detector"] = pd.Categorical(df["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df = df.sort_values("detector")

        for string_det in ["BEGe", "COAX", "ICPC", "PPC"]:
            df_string = df[df["det_type"]==string_det]
            # expo-weighted intensity
            intensity=0
            low_err=0
            upp_err=0
            for index, row in df_string.iterrows():
                intensity += row['intensity'] * row['expo']
                low_err += row['yerr_low'] * row['expo']
                upp_err += row['yerr_upp'] * row['expo']
            intensity /= df_string["expo"].sum()
            low_err /= df_string["expo"].sum()
            upp_err /= df_string["expo"].sum()

            my_label = ""
            if string_det == 1: my_label = legend_dict[cut]
            
            yerr_low_zero = (df_string["yerr_low"] == 0).all()
            yerr_upp_zero = (df_string["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                plt.errorbar(
                        str(string_det),
                        intensity,
                        yerr=[[low_err], [upp_err]],
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[my_cut], markersize=5, label=my_label
                )

            else:
                plt.errorbar(
                        str(string_det),
                        intensity,
                        yerr=[[low_err], [upp_err]],
                        marker=".", linestyle="", color=color_cut[my_cut], markersize=5, label=my_label, capsize=3
                )
            
    if peak == "K40_1461": 
        y_text = 170
    if peak == "K42_1525":
        y_text = 360
    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    #plt.xticks(rotation=90, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("exposure weighted rate [counts/(kg yr)]")
    plt.xlabel("detector type")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v],
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'],
            numpoints=1,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    if peak == "K40_1461": plt.ylim(-5,y_text+10)
    if peak == "K42_1525": plt.ylim(-5,y_text+15)
    #plt.show()
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_det_type_weighExpo.png")
        plt.savefig(f"{output_K40}/K40_det_type_weighExpo.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_det_type_weighExpo.png")
        plt.savefig(f"{output_K40}/K42_det_type_weighExpo.pdf")


################################################################################## same, but from detector type studies
RAW_path_canada = "/lfs/l1/legend/users/morella/gamma-analysis/output/vancouver-dataset/raw/"
AC_path_canada = "/lfs/l1/legend/users/morella/gamma-analysis/output/vancouver-dataset/LAr_AC/"
C_path_canada = "/lfs/l1/legend/users/morella/gamma-analysis/output/vancouver-dataset/LAr_C/"

gerda_results = "/lfs/l1/legend/users/morella/gamma-analysis/gamma-testing/attic/gerda_results/phII_phIIplus/"

gamma     = []
cut       = []
det_type  = []
intensity_canada = []
range_min_canada = []
range_max_canada = []
energy_line = []

avoid = ["e+e-_Kr85_514", "Pb212_239_Pb214_242", "Ac228_338_Pb214_352", "Ac228_965_Ac228_969", "Bi214_1120_Zn65_1125"]

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

        for my_json in json_files:
            legend = json.load(open(os.path.join(cut_path, detector_type, my_json)))
            gamma_line_name = list(legend.keys())[0]
            if gamma_line_name in avoid:
                continue
            det_type.append(detector_type)
            my_cut = ""
            if "raw" in cut_path: my_cut = "raw"
            if "LAr_C" in cut_path or "LArC" in cut_path: my_cut = "LAr_C"
            if "LAr_AC" in cut_path or "LArAC" in cut_path: my_cut = "LAr_AC"
            cut.append(my_cut)
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


for peak in peaks:
    fix, axes = plt.subplots(figsize=(10,5))
    axes.set_legend_logo(position='upper right')
    for my_cut in ["raw", "LArC", "LArAC"]:
        df_canada = single_df[(single_df.gamma == peak) & (single_df.cut == my_cut)].copy()

        df_canada["yerr_low"] = df_canada.apply(lambda row: max(row["intensity"] - row["range_min"], 0), axis=1)
        df_canada["yerr_upp"] = df_canada.apply(lambda row: max(-row["intensity"] + row["range_max"], 0), axis=1)
        df_canada.loc[df_canada["intensity"] == df_canada["yerr_low"], "yerr_low"] = 0

        # order detectors based on their position/location in the array
        valid_detectors = [detector for detector in det_list if detector in df_canada["detector"].values]
        # Convert the "detector" column to a categorical data type with the specified order
        df_canada["detector"] = pd.Categorical(df_canada["detector"], categories=valid_detectors, ordered=True)
        # Reorder the DataFrame based on the categorical order
        df_canada = df_canada.sort_values("detector")

        for string_det in ["BEGe", "COAX", "ICPC", "PPC"]:
            df_string = df_canada[df_canada["det_type"]==string_det]
            # expo-weighted intensity
            intensity=0
            low_err=0
            upp_err=0
            for index, row in df_string.iterrows():
                intensity += row['intensity'] * row['expo']
                low_err += row['yerr_low'] * row['expo']
                upp_err += row['yerr_upp'] * row['expo']
            intensity /= df_string["expo"].sum()
            low_err /= df_string["expo"].sum()
            upp_err /= df_string["expo"].sum()

            my_label = ""
            if string_det == 1: my_label = legend_dict[cut]
            
            yerr_low_zero = (df_string["yerr_low"] == 0).all()
            yerr_upp_zero = (df_string["yerr_upp"] == 0).all()

            if yerr_low_zero and yerr_upp_zero:
                plt.errorbar(
                        str(string_det),
                        intensity,
                        yerr=[[low_err], [upp_err]], 
                        #marker="_", uplims = True, linestyle="", color=color_cut[cut], markersize=5, label=legend_dict[cut], capsize=3
                        marker="x", color=color_cut[my_cut], markersize=5, label=my_label
                )

            else:
                plt.errorbar(
                        str(string_det),
                        intensity,
                        yerr=[[low_err], [upp_err]], 
                        marker=".", linestyle="", color=color_cut[my_cut], markersize=5, label=my_label, capsize=3
                )

    if peak == "K40_1461":
        y_text = 170
    if peak == "K42_1525":
        y_text = 360
    # Customize x-axis tick labels color based on the starting letter
    tick_labels = plt.xticks()[1]
    #axes.set_legend_annotation()
    #plt.xticks(rotation=90, ha='center')
    if peak == "K40_1461": plt.title(r'$^{40}$K - 1460.8 keV')
    if peak == "K42_1525": plt.title(r'$^{42}$K - 1524.7 keV')
    plt.ylabel("exposure weighted rate [counts/(kg yr)]")
    plt.xlabel("detector type")

    marker_proxy_dot = mlines.Line2D([], [], linestyle='-', marker='.', color='black', markersize=10, label='global mode +- 68%')
    marker_proxy_v = mlines.Line2D([], [], linestyle='None', marker='x', color='black', markersize=10, label='upper 90% CI limit')
    marker_proxy_dot_raw = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_v_raw = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["raw"], markersize=10, label='raw')
    marker_proxy_dot_LArC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_v_LArC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArC"], markersize=10, label='LAr C')
    marker_proxy_dot_LArAC = mlines.Line2D([], [], linestyle='-', marker='.', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    marker_proxy_v_LArAC = mlines.Line2D([], [], linestyle='None', marker='x', color=color_cut["LArAC"], markersize=10, label='LAr AC')
    axes.legend([(marker_proxy_dot_raw, marker_proxy_v_raw), (marker_proxy_dot_LArAC, marker_proxy_v_LArAC), (marker_proxy_dot_LArC, marker_proxy_v_LArC), marker_proxy_dot, marker_proxy_v],
            ["raw", 'LAr AC', "LAr C", 'global mode +- 68%', 'upper 90% limit'],
            numpoints=1,
            handler_map={tuple: HandlerTuple(ndivide=None)},
            loc = "upper left")

    if peak == "K40_1461": plt.ylim(-5,y_text+10)
    if peak == "K42_1525": plt.ylim(-5,y_text+15)
    if peak == "K40_1461":
        plt.savefig(f"{output_K40}/K40_det_type.png")
        plt.savefig(f"{output_K40}/K40_det_type.pdf")
    else:
        plt.savefig(f"{output_K42}/K42_det_type.png")
        plt.savefig(f"{output_K40}/K42_det_type.pdf")
