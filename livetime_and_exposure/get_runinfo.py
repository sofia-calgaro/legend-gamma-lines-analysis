# Get runinfo starting from skimmed hdf files.
# Pay attention to not-full runs (for some of them, some keys
# - not good for the analysis - might have been removed when 
# building the hdf file)
#
#   ! missing keys from the dataframe should be accounted
#   ! when extracting livetimes. The off time, instead, should
#   ! be accounted through meta data (and validity.jsonl time cuts).
#
# Author: Sofia Calgaro - sofia.calgaro@pd.infn.it
# Created: July 2023

import os
import json
import logging
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from matplotlib import pyplot as plt
from argparse import RawTextHelpFormatter

from legendmeta import LegendMetadata
lmeta = LegendMetadata()

import legend_data_monitor as ldm

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

def parse_json_or_dict(value):
    """Either load data stored in a JSON file or return the input dictionary."""
    try:
        # Attempt to load value as JSON
        with open(value) as json_file:
            return json.load(json_file)
    except (FileNotFoundError, json.JSONDecodeError):
        # Treat value as dictionary
        return eval(value)


def main(skimmed_path=None, version=None, prodenv=None, prodenv_version=None):
    """Get livetimes from an input file.
    
    You can run the script either from the command line, or you can use the module in another script.

    The function returns a dictionary with exposures evaluated channel by channel, run by run, period by period.
    """
    if not skimmed_path and not version and not prodenv:
        # set up the argument parser
        parser = argparse.ArgumentParser(description="Evaluate start key of runs and livetime in seconds.", formatter_class=RawTextHelpFormatter)
        
        skimmed_path_help = "Path to the hdf skimmed files."
        version_help = "Version of skimmed files."
        prodenv_help = "Path to processed lh5 files (neccessary to retrieve the first timestamp)."
        prodenv_version_help = "Version of  processed lh5 files."

        parser.add_argument("--skimmed_path", help=skimmed_path_help)
        parser.add_argument("--version", help=version_help)
        parser.add_argument("--prodenv", help=prodenv_help)
        parser.add_argument("--prodenv_version", help=prodenv_version_help)

        # parse the command-line arguments
        args = parser.parse_args()

        # get input info
        skimmed_path = args.skimmed_path
        version = args.version
        prodenv = args.prodenv
        prodenv_version = args.prodenv_version

    logger_expo.debug(f"Getting hdf files from \"{skimmed_path}\" (version: {version})")
    logger_expo.debug(f"Reading processed files from \"{prodenv}\" (version: {prodenv_version})")

    # object to fill
    run_info = {}
    run_info.update({"phy": {}, "cal": {}})

    periods = [p for p in os.listdir(os.path.join(skimmed_path, version))]

    for p in periods:
        run_info["phy"].update({p: {}})
        run_info["cal"].update({p: {}})

        runs_name = [run.split("_")[1] for run in os.listdir(os.path.join(skimmed_path, version, p)) if "hdf" in run]

        for run in runs_name:
            run_info["phy"][p].update({run: {}})
            run_info["cal"][p].update({run: {}})

            # =============================================================== CALIBRATION
            dataset = {
                "experiment": "L200", 
                "period": p,
                "type": "cal",
                "version": prodenv_version,
                "path": prodenv, 
                "runs": int(run.split("r")[-1])
            }

            _, first_timestamp, _ = ldm.utils.get_query_times(dataset=dataset)
            run_info["cal"][p][run].update({"start_key": first_timestamp})

            # =============================================================== PHYSICS
            dataset["type"] = "phy"
            _, first_timestamp, _ = ldm.utils.get_query_times(dataset=dataset)
            run_info["phy"][p][run].update({"start_key": first_timestamp})

            path = os.path.join(skimmed_path, version, p)
            data = pd.DataFrame()
            data = pd.read_hdf(f"{path}/{p}_{run}_high_level_{version}.hdf")
            data = data[data["is_pulser"] == True]
            timestamp = np.array(data.index)
            no_pulser_evts = len(np.unique(timestamp))
            #no_pulser_evts = len(data["datetime"].unique())
            livetime_in_s = no_pulser_evts * 20

            run_info["phy"][p][run].update({"livetime_in_s": livetime_in_s})

    with open(f"runinfo_{version}.json", 'w') as json_file:
        json.dump(run_info, json_file, indent=4)

    return run_info

if __name__ == "__main__":
    main()
