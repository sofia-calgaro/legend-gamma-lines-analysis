import os
import ROOT
import pandas as pd
import legend_data_monitor as ldm
from legendmeta import LegendMetadata


# specify bin width in keV
bin_width = 1 # keV
hist_range = [0, 5500]
number_bins = int((hist_range[1] - hist_range[0]) / bin_width)

# path to hdf5 data files
path = "/data1/users/krause/share/high_lvl/v02"

# list avalilable periods
periods = os.listdir(path)
print(f"Available periods: {periods}")

for p in periods:
    
    # list available runs
    runs = [run for run in os.listdir(path + "/" + p) if "hdf" in run]
    runs_name =  [run.split("_")[1] for run in os.listdir(path + "/" + p) if "hdf" in run]

    print(f"Available runs for {p}: {runs_name}")

    for r, r_name in zip(runs, runs_name):

        # load geds channelmap using legend_data_monitor
        dataset = {
            "experiment": "L200",
            "period": p,
            "type": "phy",
            "version": "",
            "path": "/data2/public/prodenv/prod-blind/tmp/auto",
            "runs": int(r_name[1:]),
        }
        geds = ldm.Subsystem("geds", dataset=dataset)
        channel_map = geds.channel_map
        
        # giving same name so that merging the two dataframes is easier
        channel_map = channel_map.rename(columns = {"channel" : "channel_id"})
        
        # initialize empty dataframe and load data
        data = pd.DataFrame()
        data = pd.read_hdf(f"{path}/{p}/{r}")

        # merge dataframes and add metadata information (detector name, string, position, etc.)
        data = data.merge(channel_map, on = "channel_id")

        # create root file in which histograms will be saved
        myfile = ROOT.TFile(f'{p}-{r_name}-spectra.root', 'RECREATE')

        print(f"Producing run {r_name}...")

        # create directories inside root files for raw, LAr AC and LAr C spectra
        raw_dir = myfile.mkdir("raw")
        lar_ac_dir = myfile.mkdir("LAr AC")
        lar_c_dir = myfile.mkdir("LAr C")
        
        # remove pulser, baseline and muon events, events not passing Quality Cuts, select events with multiplicty == 1 and load data of detectors usable for analysis
        raw = (data.is_pulser == False) & (data.is_baseline == False) & (data.is_muon_tagged == False) & (data.is_physical == True) & (data.multiplicity == 1) & (data.is_valid_channel == True)
        
        # select raw data + remove events in coincidence with LAr (LAr AC)
        lar_ac = raw & (data.is_lar_rejected == False)

        # select raw data + keep only events in coincidence with LAr (LAr C)
        lar_c = raw & (data.is_lar_rejected == True)
        
        # loop over all different kinds of data selections
        for condition, directory in zip([raw,lar_ac,lar_c],[raw_dir,lar_ac_dir, lar_c_dir]):

            directory.cd()
            
            # temporary dataframe
            tmp_data = pd.DataFrame()
            tmp_data = data[condition]
            
            # create summary histograms
            h_bege = ROOT.TH1D('BEGe', 'BEGe', number_bins, hist_range[0], hist_range[1])
            h_coax = ROOT.TH1D('COAX', 'COAX', number_bins, hist_range[0], hist_range[1])
            h_icpc = ROOT.TH1D('ICPC', 'ICPC', number_bins, hist_range[0], hist_range[1])
            h_ppc = ROOT.TH1D('PPC', 'PPC', number_bins, hist_range[0], hist_range[1])
            h_all = ROOT.TH1D('All', 'All', number_bins, hist_range[0], hist_range[1])

            # loop over HPGe detectors
            for ch in tmp_data.sort_values(by = ["location", "position"]).channel_id.unique():

                # save name, string and position
                ch_name = tmp_data[tmp_data.channel_id == ch].name.unique()[0]
                ch_string = tmp_data[tmp_data.channel_id == ch].location.unique()[0]
                ch_position = tmp_data[tmp_data.channel_id == ch].position.unique()[0]

                # create histogram
                h = ROOT.TH1D(f's{ch_string}-p{ch_position}-{ch_name}', f's{ch_string}-p{ch_position}-{ch_name}', number_bins, hist_range[0], hist_range[1])
                
                # fill single channel histogram
                for e in tmp_data[tmp_data.channel_id == ch].energy:
                    h.Fill(e)
                    if ch_name[0] == "B": h_bege.Fill(e)
                    if ch_name[0] == "C": h_coax.Fill(e)
                    if ch_name[0] == "V": h_icpc.Fill(e)
                    if ch_name[0] == "P": h_ppc.Fill(e)

                # fill all channels (full array) histogram
                for e in tmp_data.energy:
                    h_all.Fill(e)

                h.Write()
                
            h_bege.Write()
            h_coax.Write()
            h_icpc.Write()
            h_ppc.Write()
            h_all.Write()
        
        myfile.Close()