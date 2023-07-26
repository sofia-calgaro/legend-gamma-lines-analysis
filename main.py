import json
import argparse
import numpy as np
import sys
import os
import ROOT
sys.path.insert(0, './resolution')
from get_resolution import *
sys.path.insert(1, './livetime_and_exposure')
import get_exposure

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning) 

def return_config_info(config_file):
    # ...reading the config file...
    with open(config_file, "r") as file:
        config = json.load(file)
    # src path
    gamma_src_code = config["gamma-src-code"]
    # general info
    prodenv = config["dataset"]["prodenv"]
    periods = config["dataset"]["periods"]
    runs = config["dataset"]["runs"]
    detectors = config["dataset"]["detectors"]
    version = config["dataset"]["version"]
    cut = config["dataset"]["cut"]
    info = [prodenv, periods, runs, detectors, version, cut]
    # exposure info
    exposure_time_unit = config["exposure"]["time-unit"]
    run_info_path = config["exposure"]["run-info-path"]
    status = config["exposure"]["status"]
    expo = [exposure_time_unit, run_info_path, status]

    return gamma_src_code, info, expo
        
def get_histo(gamma_src_code, periods, runs, detectors, cut):
    """Combine single histograms."""
    histo_list = ROOT.std.vector('TH1D')()
    for p in periods:
        for r in runs:
            file_histo = os.path.join(gamma_src_code, "histograms/root_files", f"{p}-{r}-v02-spectra.root")
            file_root = ROOT.TFile(file_histo)
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
    return 0

def main():
    # get config file as an input
    parser = argparse.ArgumentParser(description="Main code for gamma-analysis.")
    parser.add_argument("--config", help="Path to JSON config file.")
    args = parser.parse_args()
    config_file = args.config
    print(f"You are going to inspect config={config_file}")

    # ...reading the config file...
    gamma_src_code, info, expo = return_config_info(config_file)
    print("...inspected!")
    
    #check version set by the user
    versions_avail = ["v01.05", "v01.06"]
    if not info[4] in versions_avail:
        print(f"{info[4]} is not an available version. Try among '[v01.05, v01.06]'")
        return
        
    #check cut set by the user
    cuts_avail = ['raw',  'LAr AC', 'LAr C']
    if not info[5] in cuts_avail:
        print(f"{info[5]} is not an available cut. Try among '['raw',  'LAr AC', 'LAr C']'")
        return
        
    #check periods and runs set by the user
    list_avail = []

    result_dict = {p: [i for i in range(len(info[2]) * 2)] for p in info[1]}
    exposure_det = get_exposure.main(expo[1], expo[0], str(result_dict), expo[2])
    periods_avail = list(exposure_det.keys())
    for p in periods_avail:
        runs_avail = list(exposure_det[p].keys())
        list_avail.append((p,runs_avail))

    for p in info[1]:
        if not p in np.take(list_avail,0,1):
            print(f"{p} is not an available period")
            return    
        tot_idx=len(periods_avail)
        idx=[idx for idx in range(0,tot_idx) if list_avail[idx][0]==p][0]
        for r in info[2]:
            if not r in list_avail[idx][1]:
                print(f"{p} {r} is not an available run")
                return       

            # get all the detector names
            dataset = {
                "experiment": "L200",
                "period": p,
                "type": "phy",
                "version": info[4],
                "path": info[0],
                "runs": int(r.split("r")[-1]),
            }
            geds = ldm.Subsystem("geds", dataset=dataset)
            channel_map = geds.channel_map
            all_detectors = list(channel_map.name)
        
        
    if info[3] == "single":
        for d in all_detectors:
            d=[d]
            resolution = get_resolution(config_file, d)
            histo  = get_histo(gamma_src_code, info[1], info[2], d, info[5])
            #save
            tmp_file_name = './tmp/tmp-spectra.root'
            tmp_file = ROOT.TFile(tmp_file_name, 'RECREATE')
            histo.Write()
            myfile.Close()
            #bashCommand = "make"
            #bashCommand2 = f"./runGammaAnalysis {tmp_file_name} {resolution[0]} {resolution[1]}"
            #os.system(bashCommand)
            #os.system(bashCommand2)
        return
            
            
    else:      
        resolution = get_resolution(config_file, info[3])
        histo =  get_histo(gamma_src_code, info[1], info[2], info[3], info[5])
        tmp_file_name = './tmp/tmp-spectra.root'
        tmp_file = ROOT.TFile(tmp_file_name, 'RECREATE')
        histo.Write()
        tmp_file.Close()
        bashCommand = "cd gamma-fitter-BATv100 && make"
        bashCommand2 = f"cd gamma-fitter-BATv100 && ./runGammaAnalysis {gamma_src_code} {info} {resolution[0]} {resolution[1]}"
        #os.system(bashCommand)
        #os.system(bashCommand2)
        return

    print("EOF")

if __name__=="__main__":
    main()
              
