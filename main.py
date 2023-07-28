import json
import argparse
import numpy as np
import sys
import os
import ROOT
import logging
import legend_data_monitor as ldm
sys.path.insert(0, './resolution')
from get_resolution import *
sys.path.insert(1, './livetime_and_exposure')
import get_exposure

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning) 

# -----------------------------------------------------------------------------------------
# LOGGER SETTINGS 
logger_expo = logging.getLogger(__name__)
logger_expo.setLevel(logging.DEBUG)

# create a StreamHandler to display log messages on the terminal
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)

# set the desired format for the log messages
log_formatter = logging.Formatter('%(asctime)s %(message)s')
stream_handler.setFormatter(log_formatter)

# add the StreamHandler to the logger_expo logger
logger_expo.addHandler(stream_handler)

# set the log levels for ldm.utils.logger and ldm.utils.stream_handler to shut up info from ldm code
ldm.utils.logger.setLevel(logging.WARNING) 
ldm.utils.stream_handler.setLevel(logging.WARNING)
# -----------------------------------------------------------------------------------------




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
        
def get_histo(gamma_src_code, result_dict,  detectors, cut):
    """Combine single histograms."""
    histo_list = ROOT.std.vector('TH1D')()
    for p in result_dict.keys():
        for r in result_dict[p]:
            file_histo = os.path.join(gamma_src_code, "histograms/root_files", f"{p}-{r}-v01_06-spectra.root")
            file_root = ROOT.TFile(file_histo)
            histo =  file_root.Get(f"{cut}/{detectors}")
            histo.SetDirectory(0)
            file_root.Close()
            histo_list.push_back(histo)  
    histo_tot = histo_list.at(0)
    histo_tot.SetName("Spectrum")
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
    logger_expo.debug(f"You are going to inspect config={config_file}")


    # ...reading the config file...
    gamma_src_code, output, info, expo = return_config_info(config_file)
    logger_expo.debug("...inspected!")
    
    #check version set by the user
    versions_avail = ["v01.05", "v01.06"]
    if not info[4] in versions_avail:
        logger_expo.debug(f"{info[4]} is not an available version. Try among '[v01.05, v01.06]'")
        return
        
    #check cut set by the user
    cuts_avail = ['raw',  'LAr AC', 'LAr C']
    if not info[5] in cuts_avail:
        logger_expo.debug(f"{info[5]} is not an available cut. Try among '['raw',  'LAr AC', 'LAr C']'")
        return

    #check detector set by the user
    detectors_avail = ['all',  'single', 'BEGe', 'COAX', 'ICPC', 'PPC']
    if isinstance(info[3],str) and  not info[3] in detectors_avail:
        logger_expo.debug(f"{info[3]} is not an available detector. Try among '['all',  'single', 'BEGe', 'COAX', 'ICPC', 'PPC']'")
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
            logger_expo.debug(f"{p} is not an available period")
            return
        tot_idx=len(periods_avail)
        idx=[idx for idx in range(0,tot_idx) if list_avail[idx][0]==p][0]
        run_avail= []
        for r in info[2]:
            if not r in list_avail[idx][1]:
                logger_expo.debug(f"{p} {r} is not an available run")
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
            all_detectors = list(channel_map.name)
        
        
    if info[3] == "single":
        for d in all_detectors:
            d=[d]
            output = output + "/single/" + d
            outputdir = ROOT.TNamed("outputDir",output)
            histo  = get_histo(gamma_src_code, result_dict, d, info[5])
            resolution = get_resolution(config_file, d)
            a_res = ROOT.TParameter("double")( "a_res", resolution[0] )
            b_res = ROOT.TParameter("double")( "b_res", resolution[1] )
            #store in tmp files
            tmp_directory = './tmp'
            if not os.path.exists(tmp_directory):
                os.makedirs(tmp_directory)
            tmp_file_name = f'{tmp_directory}/tmp-spectra.root'
            tmp_file = ROOT.TFile(tmp_file_name, 'RECREATE')
            outputDir.Write()
            histo.Write()
            a_res.Write()
            b_res.Write()
            tmp_file.Close()
        return
            
            
    else:
        output = output + "/" + info[3]
        outputDir = ROOT.TNamed("outputDir",output)
        histo =  get_histo(gamma_src_code, result_dict, info[3], info[5])
        resolution = get_resolution(config_file, info[3])
        a_res = ROOT.TParameter("double")( "a_res", resolution[0] )
        b_res = ROOT.TParameter("double")( "b_res", resolution[1] )
        #store resolution in json file
        resolution_file="resolution_p3p4.json"
        try:
            with open(resolution_file, 'r') as fp:
                resolution_dict = json.load(fp)
            if not info[5] in resolution_dict:
                resolution_dict[info[5]]={}
            resolution_dict[info[5]][info[3]] = resolution
            with open (resolution_file, "w") as f:
                json.dump(resolution_dict, f, indent=4)

        except FileNotFoundError:
            logger_expo.debug('Summary resolution file not found, will create a new one.')
            resolution_dict = {info[5]:{info[3]: resolution}}
            with open(resolution_file, 'w') as fp:
                json.dump(resolution_dict, fp, indent=4)
        logger_expo.debug(f'Resolution parameters are stored here: {resolution_file}')

        #store in tmp files
        tmp_directory = './tmp'
        if not os.path.exists(tmp_directory):
            os.makedirs(tmp_directory)
        tmp_file_name = f'{tmp_directory}/tmp-spectra.root'
        tmp_file = ROOT.TFile(tmp_file_name, 'RECREATE')
        outputDir.Write()
        histo.Write()
        a_res.Write()
        b_res.Write()
        tmp_file.Close()
        return

    logger_expo.debug("EOF")

if __name__=="__main__":
    main()
              
