import json
import numpy as np
import sys
import os
import ROOT
sys.path.insert(1, './resolution')
from get_resolution import *

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning) 

        
def get_histo(periods, runs, detectors, cut):
    path_histo = "/data1/users/morella/gamma-lines/legend-gamma-lines-analysis/"
    histo_list = ROOT.std.vector('TH1D')()
    for p in periods:
        for r in runs:
            file_histo = f"{p}-{r}-v02-spectra.root"
            file_root = ROOT.TFile(f"{path_histo}{file_histo}")
            histo =  file_root.Get(f"{cut}/{detectors}")
            histo.SetDirectory(0)
            file_root.Close()
            histo_list.push_back(histo)  
    histo_tot = histo_list.at(0)
    histo_tot.SetName("Sum")
    for i in range(1,histo_list.size()):
        histo_tot.Add(histo_list.at(i))
    histo_ = histo_tot.Clone()
    return histo_    

def main():
    path = '/data1/users/biancacc/gamma-lines/legend-gamma-lines-analysis/'
    config_file = f'{path}config.json'
    with open(config_file, "r") as file:
        config = json.load(file)

    periods = config["dataset"]["periods"]
    runs = config["dataset"]["runs"]
    detectors = config["dataset"]["detectors"]
    version = config["dataset"]["version"]
    cut = config["dataset"]["cut"]

    #check version set by the user
    versions_avail = ["v01.05", "v01.06"]
    if not version in versions_avail:
        print(f"{version} is not an available version. Try among '[v01.05, v01.06]'")
        return
        
    #check cut set by the user
    cuts_avail = ['raw',  'LAr AC', 'LAr C']
    if not cut in cuts_avail:
        print(f"{cut} is not an available cut. Try among '['raw',  'LAr AC', 'LAr C']'")
        return
        
    #check periods and runs set by the user
    list_avail = []
    path = '/data1/users/biancacc/gamma-lines/legend-gamma-lines-analysis/'
    file_exposure = f'{path}livetime_and_exposure/exposure_in_kg_yr_on_ac.json'
    with open(file_exposure, "r") as file:
        exposure_det = json.load(file)
    periods_avail = list(exposure_det.keys())
    for p in periods_avail:
        runs_avail = list(exposure_det[p].keys())
        list_avail.append((p,runs_avail))

    for p in periods:
        if not p in np.take(list_avail,0,1):
            print(f"{p} is not an available period")
            return    
        tot_idx=len(periods_avail)
        idx=[idx for idx in range(0,tot_idx) if list_avail[idx][0]==p][0]
        for r in runs:
            if not r in list_avail[idx][1]:
                print(f"{p} {r} is not an available run")
                return       

            # get all the detector names
            dataset = {
                "experiment": "L200",
                "period": p,
                "type": "phy",
                "version": "",
                "path": "/data2/public/prodenv/prod-blind/tmp/auto",
                "runs": int(r[1:]),
            }
            geds = ldm.Subsystem("geds", dataset=dataset)
            channel_map = geds.channel_map
            all_detectors = list(channel_map.name)
        
        
    if detectors == "single":
        for d in all_detectors:
            d=[d]
            resolution = get_resolution(version, periods, runs, d)
            histo  = get_histo(periods, runs, d, cut)
            #save
            tmp_file_name = './tmp/tmp-spectra.root'
            tmp_file = ROOT.TFile(tmp_file_name, 'RECREATE')
            histo.Write()
            myfile.Close()
            #bashCommand = "gcc runGammaAnalysis.cxx -o runGammaAnalysis"
            #bashCommand2 = f"./runGammaAnalysis {tmp_file_name} {resolution[0]} {resolution[1]}"
            #os.system(bashCommand)
            #os.system(bashCommand2)
        return
            
            
    else:      
        resolution = get_resolution(version, periods,runs, detectors)
        histo =  get_histo(periods, runs, detectors, cut)
        tmp_file_name = './tmp/tmp-spectra.root'
        tmp_file = ROOT.TFile(tmp_file_name, 'RECREATE')
        histo.Write()
        tmp_file.Close()
        #bashCommand = "gcc runGammaAnalysis.cxx -o runGammaAnalysis"
        #bashCommand2 = f"./runGammaAnalysis {tmp_file_name} {resolution[0]} {resolution[1]}"
        #os.system(bashCommand)
        #os.system(bashCommand2)
        return


if __name__=="__main__":
    main()
              