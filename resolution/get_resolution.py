import os, json, sys
import numpy as np
sys.path.insert(0, './')
import main

sys.path.insert(1, './livetime_and_exposure')
import get_exposure
from legendmeta import JsonDB
import legend_data_monitor as ldm
from legendmeta import LegendMetadata

operation ='ecal'
par ='cuspEmax_ctc_cal'
pars = 'eres_pars'

def get_resolution(config_file, detectors):
    # retrieve useful info
    gamma_src_code, output, info, expo, _ = main.return_config_info(config_file)

    file_exposure = f'exposure_in_kg_{expo[0]}'
    if isinstance(expo[2],list):
        for st in expo[2]:
            file_exposure += f"_{st}" 
    if isinstance(expo[2], str):
        file_exposure += f"_{expo[2]}" 
    file_exposure += '.json'
    with open(file_exposure, "r") as file:
        exposure_det = json.load(file)
    exposure_list = []
    resolution_list = []
    run_info = get_exposure.parse_json_or_dict(expo[1])

    for p in info[1]:
        runs_avail = exposure_det[p].keys()
        for r in info[2]:
            if not r in runs_avail:
                continue 

            # get channel map for a specific run and period
            first_timestamp = run_info["phy"][p][r]["start_key"]
            map_file = os.path.join(expo[3], expo[4], "inputs/dataprod/config")
            full_status_map = JsonDB(map_file).on(
                timestamp=first_timestamp, system="geds"
            )["analysis"]
            status_map = {key: value for key, value in full_status_map.items() if value.get('usability') in expo[2]}
            all_detectors = [det for det in status_map.keys() if "S" not in det]
            map_file = os.path.join(expo[3], expo[4], "inputs/hardware/configuration/channelmaps")
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
                path_resolutions = os.path.join(info[0], info[4], "generated/par/hit/cal", p, r)
                file= [file for file in os.listdir(path_resolutions) if file.endswith("results.json")]
                with open(os.path.join(path_resolutions, file[0]), "r") as file:
                    resolution_det = json.load(file)
                    res_det = resolution_det[c][operation][par][pars]
                resolution_list.append(res_det)


    #get the lists of the two resolution parameters
    if resolution_list != []:
        a = np.take(resolution_list,0,1)
        b = np.take(resolution_list,1,1)
    else:
        return [None, None]

    #get the average sum of the res parameters 
    par_A = weighted_average(a, exposure_list)
    par_B = weighted_average(b, exposure_list)
    weighted_pars=[par_A, par_B]

    return weighted_pars


def weighted_average(pars, weights):
    return sum(weights * pars) / sum(weights)            
