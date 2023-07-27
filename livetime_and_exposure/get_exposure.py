# Get exposure starting from a dictionary of livetimes divided run by run.
#
# Author: Sofia Calgaro - sofia.calgaro@pd.infn.it
# Created: July 2023

import os
import json
import logging
import argparse
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

def convert_time(value, unit):
    """Convert a 'value' from seconds to a given specified unit time. Choose among sec, min, hour, day, yr."""
    if unit == 'sec':
        return value
    elif unit == 'min':
        return value / 60
    elif unit == 'hour':
        return value / 3600
    elif unit == 'day':
        return value / (3600 * 24)
    elif unit == 'yr':
        return value / (3600 * 24 * 365)
    else:
        return None


def parse_json_or_dict(value):
    """Either load data stored in a JSON file or return the input dictionary."""
    try:
        # Attempt to load value as JSON
        with open(value) as json_file:
            return json.load(json_file)
    except (FileNotFoundError, json.JSONDecodeError):
        # Treat value as dictionary
        return eval(value)


def main(LT_file=None, time_unit=None, data=None, status=None):
    """Get livetimes from an input file.
    
    You can run the script either from the command line, or you can use the module in another script.

    The function returns a dictionary with exposures evaluated channel by channel, run by run, period by period.
    """
    if LT_file and time_unit and data and status:
        data = parse_json_or_dict(data)
    else:
        # set up the argument parser
        parser = argparse.ArgumentParser(description="Evaluate exposure, run by run, detector by detector, starting from input livetimes expressed in seconds.", formatter_class=RawTextHelpFormatter)
        
        livetime_help = """\
        Path to the input JSON file (or dictionary object) that stores the livetime in seconds for a given run and period and the starting key of a given run.
        In the future, this dictionary will be fully retrievable from metadata.
        The content of the dictionary has to be of the following type:

        {
        "phy": {
            "p03": {
                "r000": {
                    "start_key": "20230312T043356Z",
                    "livetime_in_s": 443740
                },
                "r001": {
                    "start_key": "20230318T015140Z",
                    "livetime_in_s": 550760
                },
                "r002": {
                    "start_key": "20230324T205907Z",
                    "livetime_in_s": 554660
                },
                ...
            },
            ...
        }

        Note: If you use a dictionary from the command line, be sure to enclose the dictionary argument within single quotes (') 
        to ensure that the entire dictionary is treated as a single string.
        """

        time_unit_help = """\
        String (among 'sec', 'min', 'hour', 'day', 'yr') to get the exposure in a given time format. Mass is always computed in kg in the exposure's evaluation.
        """

        data_help = """\
        Path to the input JSON file (or dictionary object) that stores periods and runs for which we want to compute the exposure.
        The content of the dictionary has to be of the following type:

        {
            "p03": [0, 1, 2],
            "p04": [0, 1]
        }

        Note: If you use a dictionary from the command line, be sure to enclose the dictionary argument within single quotes (') 
        to ensure that the entire dictionary is treated as a single string.
        """

        status_help = """\
        Statuses (among 'on', 'off', 'ac') to include when evaluating the exposure.
        The argument has to be of the following type:
        
        1. on = to include only ON detectors 
        2. on off = to include only ON+OFF detectors
        3. on off ac = to include only ON+OFF+AC detectors
        """

        parser.add_argument("--livetime", help=livetime_help)
        parser.add_argument("--time_unit", help=time_unit_help)
        parser.add_argument("--data", help=data_help)
        parser.add_argument("--status", nargs="+", help=status_help)

        # parse the command-line arguments
        args = parser.parse_args()

        # get input info
        LT_file = args.livetime
        time_unit = args.time_unit
        data = parse_json_or_dict(args.data)
        status = args.status

    logger_expo.debug(f"Getting livetime values from \"{LT_file}\"")
    logger_expo.debug(f"You are going to evaluate exposure in kg*{time_unit}")
    logger_expo.debug(f"You are going to inspect following periods and runs: {data}")
    logger_expo.debug(f"You are going to inspect detectors that have the following status: {status}")
    logger_expo.debug(f"Individual exposures will be stored here: \"exposure_in_kg_{time_unit}.json\"")
    logger_expo.debug(f"Summary exposures (and masses) will be stored here: \"output.log\"")

    # get hardware map (for masses)
    dets_map = lmeta.hardware.detectors.germanium.diodes

    # Remove the file if it already exists
    if os.path.exists('output.log'):
        os.remove('output.log')

    # get livetimes
    run_info = parse_json_or_dict(LT_file)

    expo_all_periods_runs = {}

    # inspect each periods individually
    for period in data.keys():
        logger_expo.info(f"... inspecting period {period}")
        # available runs
        runs = data[period]
        expo_all_periods_runs.update({period: {}})

        # inspect each run individually
        for run in runs:
            r_run = f"{str(run).zfill(3)}"
            logger_expo.info(f"...... {r_run}")
            expo_all_periods_runs[period].update({r_run: {}})

            # get channel map for a specific run and period
            first_timestamp = run_info["phy"][period][r_run]["start_key"]
            full_map = lmeta.dataprod.config.on(
                timestamp=first_timestamp, system="geds"
            )["analysis"]

            # keep detectors of for a given specified status
            channel_map = {key: value for key, value in full_map.items() if value.get('usability') in status}

            # get a list of detectors (exclude SiPMs)
            det_names = [det for det in channel_map.keys() if "S" not in det]

            # ===================================================
            # Get livetime from input dictionary
            # ===================================================
            livetime_run = run_info["phy"][period][r_run]["livetime_in_s"]

            tot_expo = 0
            tot_expo_coax = 0
            tot_expo_bege = 0
            tot_expo_icpc = 0
            tot_expo_ppc = 0
            tot_mass_coax = 0
            tot_mass_bege = 0
            tot_mass_icpc = 0
            tot_mass_ppc = 0

            # ===================================================
            # Get exposure - channel by channel
            # ===================================================
            for det_name in det_names:
                mass_in_kg = dets_map[det_name]["production"]["mass_in_g"] / 1000
                
                ch_expo = mass_in_kg * convert_time(livetime_run, time_unit)
                # stored exposure - for later output
                tot_expo += ch_expo
                tot_expo_coax += ch_expo if "C" == det_name[0] else 0
                tot_expo_bege += ch_expo if "B" == det_name[0] else 0
                tot_expo_icpc += ch_expo if "V" == det_name[0] else 0
                tot_expo_ppc += ch_expo if "P" == det_name[0] else 0
                # stored mass - for later output
                tot_mass_coax += mass_in_kg if "C" == det_name[0] else 0
                tot_mass_bege += mass_in_kg if "B" == det_name[0] else 0
                tot_mass_icpc += mass_in_kg if "V" == det_name[0] else 0
                tot_mass_ppc += mass_in_kg if "P" == det_name[0] else 0

                # save exposure for a given channel
                expo_all_periods_runs[period][r_run].update({det_name: ch_expo})
            
            # save summary exposure in output file
            logging.basicConfig(filename='output.log', level=logging.INFO, format='%(message)s')
            logging.info(f"\nTotal exposure for {period}-{run} is: {round(tot_expo,2)} kg*{time_unit}")
            logging.info(f"--- Total exposure for COAX: {round(tot_expo_coax,2)} kg*{time_unit} (m_tot={round(tot_mass_coax,2)} kg)")
            logging.info(f"--- Total exposure for BEGe: {round(tot_expo_bege,2)} kg*{time_unit} (m_tot={round(tot_mass_bege,2)} kg)")
            logging.info(f"--- Total exposure for ICPC: {round(tot_expo_icpc,2)} kg*{time_unit} (m_tot={round(tot_mass_icpc,2)} kg)")
            logging.info(f"--- Total exposure for PPC: {round(tot_expo_ppc,2)} kg*{time_unit} (m_tot={round(tot_mass_ppc,2)} kg)")

    # save the single-channel exposure dictionary to a JSON file
    output_json_file = f'exposure_in_kg_{time_unit}'
    if isinstance(status,list):
        for st in status:
            output_json_file += f"_{st}" 
    if isinstance(status, str):
        output_json_file += f"_{status}" 
    output_json_file += '.json'
    with open(output_json_file, 'w') as json_file:
        json.dump(expo_all_periods_runs, json_file, indent=4)

    return expo_all_periods_runs

if __name__ == "__main__":
    main()
