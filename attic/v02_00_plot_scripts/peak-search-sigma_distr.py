#from legend_plot_style import LEGENDPlotStyle as lps

import os
import json
import re
import numpy as np
import pandas as pd

from legendmeta import JsonDB
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

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
                    "e+e-_511"    :    511.0,
                    "Kr85_514"    :    514.0,
                    "Pb212_239"   :    238.6,
                    "Pb214_242"   :    242.0,
                    "Ac228_338"   :    338.3,
                    "Pb214_352"   :    351.9,
                    "Ac228_965"   :    964.8,
                    "Ac228_969"   :    969.0,
                    "Bi214_1120"  :    1120.3,
                    "Zn65_1125"   :    1125.0,
                    "Ac228_1588"  :    1588 }

E_gamma = [238.6, 295.224, 338.0, 352.0, 434.0, 478.3, 514.0, 583.2, 609.3, 662.0, 727.0, 911.2, 969.0, 1001.3, 1173.2, 1238.0, 1332.5, 1408.0, 1460.8, 1524.6]

dataset_type = "All"
raw_all_folder = f"../output_gamma_rates/vancouver_silver_peak_search_LArAC/{dataset_type}/"
output_folder = "peak_search"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

import ROOT

energies = []
s_68 = []
global_modes = []
rng_min = []
rng_max = []
gm_norm_to_sigma68 = []
gm_diff_to_sigma68 = []
p0_constant = []
p1_slope = []
excesses_3sigma = []

all_json_files = os.listdir(raw_all_folder)
all_energies = [(f.split('.json')[0]).split('histo.')[-1] for f in all_json_files if 'json' in f]
all_energies_float = sorted([float(en) for en in all_energies])
all_energies = [str(en)  for en in all_energies_float]
print(f"You'll inspect energies from {all_energies[0]} to {all_energies[-1]}")

for E0 in all_energies:
    json_file = os.path.join(raw_all_folder, f"histo.{E0}.json")
    with open(json_file, 'r') as file:
        my_json = json.load(file)

    if "mode" in my_json[E0]["fit_parameters"]["line"]["intensity0_in_cts"].keys():
        #print(E0)
        fwhm = my_json[E0]["fit_parameters"]["line"]["fwhm0_in_keV"]["value"]
        gm = my_json[E0]["fit_parameters"]["line"]["intensity0_in_cts"]["mode"] 
        E0_counts_L68 = my_json[E0]["fit_parameters"]["line"]["intensity0_in_cts"]["range_min"]
        E0_counts_U68 = my_json[E0]["fit_parameters"]["line"]["intensity0_in_cts"]["range_max"] 
        global_modes.append(gm)
        rng_min.append(E0_counts_L68)
        rng_max.append(E0_counts_U68)
        p0 = my_json[E0]["fit_parameters"]["background"]["par0.0_in_cts"]["value"]
        p1 = my_json[E0]["fit_parameters"]["background"]["par1_in_cts/keV"]["value"]
        p0_constant.append(p0)
        p1_slope.append(p1)

        root_file = f"histo_marginalized.{E0}.root"
        file_root = ROOT.TFile.Open(os.path.join(raw_all_folder, root_file), "READ")
        if not os.path.join(raw_all_folder, root_file):
            print("error", E0)
            continue
        h = file_root.Get("h1_histo_fitter_model_parameter_intensity0")
        h_clone = file_root.Get("h1_histo_fitter_model_parameter_intensity0")

        binmax = h.GetMaximumBin()
        GlobalMode = gm
        for k in range(0,h.GetNbinsX(),1):
            h.SetBinContent(k, h.GetBinContent(k+h.FindBin(GlobalMode)-1)) # shift to left        
        p = [0, 0.3415, 0.46, 0.50, 0.84, 0.90, 0.95]
        q = np.zeros(len(p), dtype=float)
        h_clone.GetQuantiles(len(p), q, np.array(p, dtype=float)) # GetQuantiles(n_division, quantiles, probsum[%])
        qt_34 = q[1] - q[0]

        # delete histos and close file
        h = None
        h_clone = None
        file_root.Close()

        if GlobalMode/qt_34 < 1:
            sigma_68 = qt_34
        else:
            sigma_68 = ( E0_counts_U68 - E0_counts_L68 )*0.5

        if GlobalMode/sigma_68 > 3:
            print(f"Excesss for {E0} at {GlobalMode/sigma_68}")
            excesses_3sigma.append(float(E0))

        energies.append(float(E0))
        s_68.append(sigma_68)
        gm_norm_to_sigma68.append(GlobalMode/sigma_68)
        gm_diff_to_sigma68.append(GlobalMode - 5 * sigma_68)

data = [
    [energies, s_68, 5 * s_68, global_modes, gm_diff_to_sigma68, gm_norm_to_sigma68],
]
import pandas

plt.figure(1, figsize=(12,6))
identified_lines = [238.6, 295.224, 351.932, 478.4, 513.997, 583.187, 661.657, 911.196]
for en in gamma_energy.values():
    if en > float(E0): continue
    plt.axvline(x=en, color="grey", linestyle="--", linewidth=1)
    if en == 238.6: plt.axvline(x=en, color="grey", linestyle="--", linewidth=1, label="Already inspected $\gamma$")
#for en in gamma_energy.values():
#    if float(en) > np.array(energies).max():
#        continue
#    plt.axvline(x=float(en), color='gray', linestyle='--')
plt.axhline(y = 3, color="darkorange", label=r'$3 \sigma$ threshold', linewidth=1)	
#plt.axhline(y = 5, color="firebrick", label=r'$5 \sigma$ thr.', linewidth=1)	

plt.plot(energies, gm_norm_to_sigma68, '-', color="royalblue", linewidth=1)	
plt.xlabel(r'Energy [keV]')
plt.ylabel(r'Significance, $\sigma$')
plt.legend(ncol=1, loc = "upper right", frameon=True, facecolor='white', framealpha=1) 
plt.title(f"Peak search over {dataset_type} spectrum")
plt.tight_layout()
plt.savefig(f"{output_folder}/energy_excesses.pdf")
plt.savefig(f"{output_folder}/energy_excesses.png")
plt.ylim(-1,10)
plt.savefig(f"{output_folder}/energy_excesses_zoom.pdf")
plt.savefig(f"{output_folder}/energy_excesses_zoom.png")
plt.close()

plt.figure(2, figsize=(12,6))
plt.plot(energies, p0_constant, color='red')
plt.xlabel(r'Energy [keV]')
plt.ylabel('p0, constant [cts]')
plt.title(f"Peak search over {dataset_type} spectrum")
plt.tight_layout()
plt.savefig(f"{output_folder}/energies_vs_p0.pdf")
plt.savefig(f"{output_folder}/energies_vs_p0.png")
plt.close()


plt.figure(3, figsize=(12,6))
plt.plot(energies, p1_slope)
plt.xlabel(r'Energy [keV]')
plt.ylabel('p1, slope [cts/keV]')
plt.title(f"Peak search over {dataset_type} spectrum")
plt.tight_layout()
plt.savefig(f"{output_folder}/energies_vs_p1.pdf")
plt.savefig(f"{output_folder}/energies_vs_p1.png")
plt.close()

