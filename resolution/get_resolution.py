import os, json, sys
import numpy as np
sys.path.insert(1, './')
import main

import legend_data_monitor as ldm
from legendmeta import LegendMetadata

operation ='ecal'
par ='cuspEmax_ctc_cal'
pars = 'eres_pars'

def get_resolution(config_file, detectors):
    # retrieve useful info
    gamma_src_code, info, expo = main.return_config_info(config_file)

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

    for p in info[1]:
        runs_avail = exposure_det[p].keys()
        for r in info[2]:
            if not r in runs_avail:
                continue 

            # load geds channelmap using legend_data_monitor
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
            if detectors == "all":
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
                detector_list = detectors

            for d in detector_list:
                if not d in all_detectors: #for list of detectors set by the user
                    continue
                if d == "B00089D":  #set to on but it is off
                    continue

                # get exposure
                exp_det = exposure_det[p][r][d]
                exposure_list.append(exp_det)

                #get resolution
                c=channel_map.channel[channel_map.name==d]
                c="ch"+str(c.values[0])
                path_resolutions = os.path.join(info[0], info[4], "generated/par/hit/cal", p, r)
                file= [file for file in os.listdir(path_resolutions) if file.endswith("results.json")]
                with open(os.path.join(path_resolutions, file[0]), "r") as file:
                    resolution_det = json.load(file)
                    res_det = resolution_det[c][operation][par][pars]
                resolution_list.append(res_det)

    #get the lists of the two resolution parameters
    a = np.take(resolution_list,0,1)
    b = np.take(resolution_list,1,1)

    #get the average sum of the res parameters 
    par_A = weighted_average(a, exposure_list)
    par_B = weighted_average(b, exposure_list)
    weighted_pars=[par_A, par_B]

    return weighted_pars


def weighted_average(pars, weights):
    return sum(weights * pars) / sum(weights)            