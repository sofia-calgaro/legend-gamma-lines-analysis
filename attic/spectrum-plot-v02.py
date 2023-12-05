from legend_plot_style import LEGENDPlotStyle as lps
import os
import re
import json
import ROOT
import numpy as np
import pandas as pd
import ROOT
import argparse
import pandas as pd
import lgdo.lh5_store as lh5
import legend_data_monitor as ldm
from legendmeta import JsonDB
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
#lps.use('legend_talks')

output_folder = "spectrum_plots"

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

def return_final_histo(my_periods, my_runs, my_root_files_path, cut, det):
    histo_list = ROOT.std.vector('TH1D')()
    for idx,period in enumerate(my_periods):
        for run in my_runs[idx]:
            
            file_root = ROOT.TFile(os.path.join(my_root_files_path, f"{period}-{run}-v02.00-spectra.root"))
            histo_flag = f"{cut}/{det}"
            print(histo_flag)
            histo =  file_root.Get(histo_flag)
            if not histo:
               print(f"no file for {period}-{run}")
               file_root.Close()
               continue
            histo.SetDirectory(0)
            file_root.Close()
            histo_list.push_back(histo)

    if len(histo_list) == 0:
        print(my_periods, my_runs, my_root_files_path, cut, det)
        return None
    histo_tot = histo_list.at(0)
    histo_tot.SetName("Spectrum")
    for i in range(1,histo_list.size()):
        histo_tot.Add(histo_list.at(i))
    final_histo = histo_tot.Clone()

    return final_histo


def return_histo_2(my_periods, my_runs, my_root_files_path, cut, det):
    final_histo  = ROOT.TH1D("Spectrum", "", 5000,0,5000)
    for idx,period in enumerate(my_periods):
        for run in my_runs[idx]:

            file_root = ROOT.TFile(os.path.join(my_root_files_path, f"{period}-{run}-v02.00-spectra.root"))
            histo_flag = f"{cut}/{det}"
            print(histo_flag)
            histo =  file_root.Get(histo_flag)
            if not histo:
               print(f"no file for {period}-{run}")
               file_root.Close()
               continue
            final_histo.Add(histo)

    return final_histo


full_path = "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/src/root_files"

taup_root_files_path = os.path.join(full_path, "p03_p07_v02_00_1keVbinning")
taup_periods = ["p03", "p04"]
taup_runs = [["r000", "r001", "r002", "r003", "r004", "r005"], ["r000","r001", "r002", "r003"]]
taup_dict_exposure = {"ICPC": 9.35, "BEGe": 2.1, "COAX": 1.86, "PPC": 1.42}
taup_dict_exposure.update({"All": sum(taup_dict_exposure.values())})

canada_root_files_path = os.path.join(full_path, "p03_p07_v02_00_1keVbinning")
canada_periods = ["p03", "p04", "p06", "p07"]
canada_runs = [["r000", "r001", "r002", "r003", "r004", "r005"], ["r000","r001", "r002", "r003"], ["r000", "r001", "r002", "r003", "r004", "r005"], ["r001", "r002", "r003", "r004", "r005", "r006", "r007"]]
canada_dict_exposure = {"ICPC": 27.2, "BEGe": 6.3, "COAX": 5.34, "PPC": 4.9}
canada_dict_exposure.update({"All": sum(canada_dict_exposure.values())})

legend_color_cut = {"raw" : "dodgerblue", "LArAC": "nvay", "LArC" : "orange"}
label_cut = {"raw": "before LAr", "LAr AC": "after LAr", "LAr C" : "with LAr"}

bin_width = 2 # keV
n_bins = int(5000 / bin_width)

# ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL ALL
detectors = ["All", "ICPC", "COAX", "BEGe", "PPC"]
cuts = ["raw", "LAr C", "LAr AC"]
my_cut = {"raw": "raw", "LAr C": "LArC", "LAr AC": "LArAC"}
my_x_zooms = [[200,5000], [200, 800], [800,1600], [1600,2800]]
y_txt = [100, 100, 100, 3]
list_custom_lines = [k for k in gamma_isotope.keys()]
list_custom_lines = sorted(list_custom_lines, key=lambda label: gamma_energy[label])

for idx_zoom,x_zoom in enumerate(my_x_zooms):
    for det in detectors:
        for cut in cuts:
            plt.figure(figsize = (12,6))

            # TAUP 
            taup_histo = return_final_histo(taup_periods, taup_runs, taup_root_files_path, cut, det)
            if not taup_histo: continue
            taup_histo.Rebin(bin_width)
            taup_no_bins = taup_histo.GetNbinsX()
            taup_bin_cntnt = []
            for i in range(0, taup_no_bins):
                taup_bin_cntnt.append(taup_histo.GetBinContent(i))
            taup_bin_cntnt = np.array(taup_bin_cntnt)/taup_dict_exposure[det]
            plt.hist(np.arange(-1, taup_no_bins-1, 1)*bin_width, weights = taup_bin_cntnt, bins = taup_no_bins, histtype = "step", linewidth = 1, color = "dodgerblue", label = f"{det} - TAUP")
            #my_file = ROOT.TFile(f"{output_folder}/{my_cut[cut]}-{det}-spectrum-{x_zoom[0]}_{x_zoom[1]}.ROOT", "RECREATE")
            #taup_histo.Write()
            #my_file.Close()

            # Vancouver canada_
            canada_histo = return_final_histo(canada_periods, canada_runs, canada_root_files_path, cut, det)
            if not canada_histo: continue
            canada_histo.Rebin(bin_width)
            canada_no_bins = canada_histo.GetNbinsX()
            canada_bin_cntnt = []
            for i in range(0, canada_no_bins):
                canada_bin_cntnt.append(canada_histo.GetBinContent(i))
            canada_bin_cntnt = np.array(canada_bin_cntnt)/canada_dict_exposure[det]
            plt.hist(np.arange(-1, canada_no_bins-1, 1)*bin_width, weights = canada_bin_cntnt, bins = canada_no_bins, histtype = "step", linewidth = 1, color = "orange", label = f"{det} - Vancouver")

            plt.legend(loc = "lower left")
            plt.xlabel("energy [keV]")
            plt.ylabel(f"counts / {bin_width} keV / (kg yr)")
            plt.yscale("log")
            plt.xlim(x_zoom)
            plt.ylim([0.01, 2e2])
            plt.title(f"Spectrum {label_cut[cut]}")
            plt.tight_layout()

            print("taup:", taup_dict_exposure[det], "\tvancouver:", canada_dict_exposure[det])

            if idx_zoom>0:
                for k_idx,k in enumerate(list_custom_lines):
                    if gamma_energy[k] < x_zoom[0] or gamma_energy[k] > x_zoom[1]: continue
                    if k == "e+e-_511": continue

                    plt.axvline(x=gamma_energy[k], color='gray', linestyle="--")
                    scarto = 0
                    if idx_zoom != 10:
                        if k_idx % 2 == 0:
                            if y_txt[idx_zoom]==100:
                                scarto = 30
                            else:
                                scarto = 1
                    if k == "e+e-": scarto = 0
                    dx = 0.5
                    if gamma_energy[k] in [238.6, 511, 964.8, 1120.3, 295.2]:
                        dx = -30
                    plt.text(gamma_energy[k]+dx, y_txt[idx_zoom]-scarto, gamma_isotope[k], fontsize=14)

            if not os.path.exists(output_folder):
               os.makedirs(output_folder)

            zoom_str = f"{x_zoom[0]}_{x_zoom[1]}"

            plt.savefig(f"{output_folder}/{my_cut[cut]}-{det}-spectrum-{zoom_str}.pdf")
            plt.savefig(f"{output_folder}/{my_cut[cut]}-{det}-spectrum-{zoom_str}.png")
