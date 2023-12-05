#from legend_plot_style import LEGENDPlotStyle as lps
import os
import json
import re
import numpy as np
import pandas as pd
import numpy as np
from legendmeta import JsonDB
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from legendmeta import LegendMetadata
import json
import pygama.math.histogram as pgh
import pygama.math.peak_fitting as pgf
from pygama.pargen.energy_optimisation import fwhm_slope
from legendmeta import JsonDB
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

gamma_energy = {    "K42_1525"  :    1524.7,
                    "K40_1461"  :    1460.8 }

periods = ["p03"]#, "p04", "p06", "p07"]
runs = [["r000"]]#, "r001", "r002", "r003", "r004", "r005"], ["r000","r001", "r002", "r003"], ["r000", "r001", "r002", "r003", "r004", "r005"], ["r001", "r002", "r003", "r004", "r005", "r006", "r007"]]  

for idx,period in enumerate(periods):
    for run in runs[idx]:
        print(f"... {period} -- {run}")
        det_dict = {}
        for peak in gamma_energy.keys():
            for idx,period in enumerate(periods):
                for run in runs[idx]:
                    raw_folder = f"../output/{period}_{run}/single"
                    detectors = os.listdir(raw_folder)

                    for det in detectors:
                        json_file = os.path.join(raw_folder, det, f"histo.{peak}.json")

                        with open(json_file, 'r') as file:
                            my_json = json.load(file)

                        fwhm = my_json[peak]["fit_parameters"]["line"]["fwhm0_in_keV"]["value"]
                        fwhm_err = my_json[peak]["fit_parameters"]["line"]["fwhm0_in_keV"]["err"]
                        centroid = my_json[peak]["fit_parameters"]["line"]["mean0_in_keV"]["value"] 
                        centroid_err = my_json[peak]["fit_parameters"]["line"]["mean0_in_keV"]["err"]
                        line = my_json[peak]["fit_parameters"]["line"]["mean0_in_keV"]["line"]
                        residual = centroid - line
                        residual_err = centroid_err
                        if det not in det_dict.keys():
                            det_dict.update({det: {peak: {"fwhm": fwhm, "fwhm_err": fwhm_err, "centroid": centroid, "centroid_err": centroid_err, "residual": residual, "residual_err": residual_err}}})
                        else:
                            det_dict[det].update({peak: {"fwhm": fwhm, "fwhm_err": fwhm_err, "centroid": centroid, "centroid_err": centroid_err, "residual": residual, "residual_err": residual_err}})


        fig0, axis0 = plt.subplots(figsize=(15,6))
        fig1, axis1 = plt.subplots()
        dx = 0.1
        plot_range = (-2,2)
        fit_range = (-1,1)
        cols = ['limegreen', 'dodgerblue', 'darkgreen', 'navy']

        map_file = os.path.join("/lfs/l1/legend/data/public/prodenv/prod-blind/ref/", "v02.00", "inputs/hardware/configuration/channelmaps")
        first_key = '20230311T235840Z' # start of p03
        chmap  = JsonDB(map_file).on(timestamp=first_key)
        channels = [chn for chn in chmap.keys() if chmap[chn]['system']=='geds'] # select HPGe channels
        my_channels = [ch for ch in channels if ch in det_dict.keys()]
        strings = np.array([chmap[chn]['location']['string'] for chn in my_channels])
        string_lines = [k for k, s in enumerate(strings[:-1]) if s < strings[k+1]]

        for j, peak in enumerate(gamma_energy.keys()):
                residuals = np.array([det_dict[chn][peak]['residual'] for chn in det_dict.keys()])
                residuals_err = np.array([det_dict[chn][peak]['residual_err'] for chn in det_dict.keys()])
                axis0.errorbar(my_channels, residuals, fmt='o',color=cols[j],ms=4,ls='',label=f'{peak} keV')
                axis0.fill_between(my_channels, residuals - residuals_err, residuals + residuals_err, color=cols[j], alpha=0.5)
        
                for line in string_lines: axis0.axvline(line+0.5,color='k',ls=':')
                # residual distribution
                hist, bins, var = pgh.get_hist(residuals, range=plot_range, dx=dx)
                bin_cs = (bins[:-1] + bins[1:]) / 2
                guess = (0,0.1,hist.max())
                #bounds = (None,(0,None),(0,None))
                #pars, errs = ut.fit_binned(pgf.gauss_amp, hist, bins, var, fit_range, cost_func = 'Least Squares', guess = guess, bounds = bounds)
                #fwhm, fwhm_err = pars[1]*2.355, errs[1]*2.355
                axis1.plot(bin_cs,hist,ds='steps',c=cols[j])
                axis1.fill_between(bin_cs, hist, step="pre", color=cols[j], alpha=0.5, label=f'{peak} keV')# keV\nmean = {pars[0]:.2f} ± {errs[0]:.2f} keV\nFWHM = {fwhm:.2f} ± {fwhm_err:.2f} keV')
                #axis0.set_legend_logo(position='upper left',scaling_factor=8)
                #axis0.set_legend_annotation(month)
                axis0.legend(loc='upper right')

        plt.sca(axis0)
        plt.xticks(rotation=90,fontsize=12)
        plt.sca(axis1)
        axis0.set_ylabel('$E_{cal} - E_{true}$ (keV)')
        axis1.set_xlabel('$E_{cal} - E_{true}$ (keV)')
        axis1.set_ylabel(f'Counts/({dx:.1f} keV)')
        #axis1.set_legend_logo(position='upper left')
        #axis1.set_legend_annotation(month)
        axis1.legend(loc='upper right', labelspacing = 1.0)
        #axis1.set_ylim(0,12)
        folder_path = f"{period}_{run}"
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        fig0.savefig(os.path.join(folder_path, "phy_residuals.png"),dpi=300, bbox_inches='tight')
        fig1.savefig(os.path.join(folder_path, 'phy_residuals_hist.png'),dpi=300, bbox_inches='tight')

        with open(os.path.join(folder_path, "residuals.json"), 'w') as json_file:
                json.dump(det_dict, json_file, indent=4)
