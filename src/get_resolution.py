import os, json, sys
import numpy as np
sys.path.insert(0, './')
import main
from legendmeta import JsonDB

from legendmeta import LegendMetadata
lmeta = LegendMetadata()

par ='cuspEmax_ctc_cal'
# pars = 'eres_pars'
pars = 'eres_linear'

def get_resolution(config_file, detectors):
    # retrieve useful info
    _, _, info, status, _ = main.return_config_info(config_file)

    file_exposure = f'src/settings/exposure_in_kg_yr'
    if isinstance(status,list):
        for st in status:
            file_exposure += f"_{st}" 
    if isinstance(status, str):
        file_exposure += f"_{status}" 
    file_exposure += '.json'
    with open(file_exposure, "r") as file:
        exposure_det = json.load(file)
    exposure_list = []
    resolution_list = []
    run_info = lmeta.dataprod.runinfo

    for idx_p,p in enumerate(info[1]):
        runs_avail = exposure_det[p].keys()
        for r in info[2][idx_p]:
            if not r in runs_avail:
                continue 

            # get channel map for a specific run and period
            first_timestamp = run_info[p][r]["phy"]["start_key"]
            map_file = os.path.join(info[0], info[4], "inputs/dataprod/config")
            full_status_map = JsonDB(map_file).on(
                timestamp=first_timestamp, system="geds"
            )["analysis"]
            status_map = {key: value for key, value in full_status_map.items() if value.get('usability') in status}
            all_detectors = [det for det in status_map.keys() if "S" not in det]
            map_file = os.path.join(info[0], info[4], "inputs/hardware/configuration/channelmaps")
            channel_map = JsonDB(map_file).on(timestamp=first_timestamp)

            if detectors == "All":
                detector_list = all_detectors
            elif detectors == "BEGe":
                detector_list= [idx for idx in all_detectors if idx[0] == "B"]
            elif detectors == "ICPC":
                detector_list= [idx for idx in all_detectors if idx[0] == "V"]
            elif detectors == "COAX":
                detector_list= [idx for idx in all_detectors if idx[0] == "C"]
            elif detectors == "PPC":
                detector_list = [idx for idx in all_detectors if idx[0] == "P"]
            else:
                detector_list = [detectors]

            for d in detector_list:
                if len(detector_list) == 1 and isinstance(detector_list, list):
                    d = d[0] if isinstance(d, list) else d
                print(".....", d)
                if d not in all_detectors: # for list of detectors set by the user
                    continue
                if d not in list(exposure_det[p][r].keys()): # skip a detector if not present in the exposure list
                    continue

                # get exposure
                exp_det = exposure_det[p][r][d]
                exposure_list.append(exp_det)

                #get resolution
                c="ch"+str(channel_map[d]['daq']['rawid'])
                #path_resolutions = os.path.join(info[0], info[4], "generated/par/hit/cal", p, r)
                path_resolutions = os.path.join(info[0], info[4], "generated/par/pht/cal", p, r)
                file = [file for file in os.listdir(path_resolutions) if file.endswith(".json")]

                with open(os.path.join(path_resolutions, file[0]), "r") as file:
                    resolution_det = json.load(file)
                    if info[4]=="v02.00" and p=='p07' and d in ["V07646A", "V05261A"]:
                        # otherwise, resolution for these detectors is NaN in v02.00 partitions 
                        res_det = resolution_det[c]["results"]['ecal'][par][pars]["parameters"]
                    else:
                        res_det = resolution_det[c]["results"]['partition_ecal'][par][pars]["parameters"]
                
                if not np.isnan(res_det['a']) and not np.isnan(res_det['b']):
                    resolution_list.append(res_det)


    #get the lists of the two resolution parameters
    if resolution_list != []:
        # a = np.take(resolution_list,0,1)
        # b = np.take(resolution_list,1,1)
        a = np.array([float(dictionary["a"]) for dictionary in resolution_list])
        b = np.array([float(dictionary["b"]) for dictionary in resolution_list])
    else:
        return [None, None]

    #get the average sum of the res parameters 
    par_A = weighted_average(a, exposure_list)
    par_B = weighted_average(b, exposure_list)
    weighted_pars=[par_A, par_B]

    return weighted_pars


def weighted_average(pars, weights):
    return sum(weights * pars) / sum(weights)            
