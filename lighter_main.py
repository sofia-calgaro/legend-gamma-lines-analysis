import json
import argparse
import numpy as np
import sys
import os
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
    # output path
    output = config["output"]
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

    return gamma_src_code, output, info, expo
     

def main():
    # get config file as an input
    parser = argparse.ArgumentParser(description="Main code for gamma-analysis.")
    parser.add_argument("--config", help="Path to JSON config file.")
    args = parser.parse_args()
    config_file = args.config
    print(f"You are going to inspect config={config_file}")


    # ...reading the config file...
    gamma_src_code, output, info, expo = return_config_info(config_file)
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

    with open(expo[1], "r") as file:
        run_info = json.load(file)
    periods_avail = list(run_info["phy"].keys())
    for p in periods_avail:
        runs_avail = list(run_info["phy"][p].keys())
        list_avail.append((p,runs_avail))

    result_dict = {}

    for p in info[1]:
        if not p in np.take(list_avail,0,1):
            print(f"{p} is not an available period")
            return
        tot_idx=len(periods_avail)
        idx=[idx for idx in range(0,tot_idx) if list_avail[idx][0]==p][0]
        run_avail= []
        for r in info[2]:
            if not r in list_avail[idx][1]:
                print(f"{p} {r} is not an available run")
                #return
            else:
                run_avail.append(r)
        result_dict[p] = run_avail   

    exposure_det = get_exposure.main(expo[1], expo[0], str(result_dict), expo[2])
    periods_avail = list(exposure_det.keys())

    for p in result_dict.keys():
        for r in result_dict[p]:
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
            channel_map = channel_map[channel_map.status != "off"]
            all_detectors = list(channel_map.name)
        
    total_resolution = {}
        
    if info[3] == "single":
        for d in all_detectors:
            print(".....", d)
            total_resolution.update({d: {}})
            d=[d]
            resolution = get_resolution(config_file, d)
            total_resolution[d[0]].update({"a_res": resolution[0], "b_res": resolution[1]})
            
            
    else:
        total_resolution.update({info[3]: {}})
        resolution = get_resolution(config_file, info[3])
        total_resolution[info[3]].update({"a_res": resolution[0], "b_res": resolution[1]})

    with open("weighted_resolutions.json", 'w') as json_file:
        json.dump(total_resolution, json_file, indent=4)

    print("EOF")
    return

if __name__=="__main__":
    main()
