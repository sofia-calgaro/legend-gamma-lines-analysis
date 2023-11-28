import os
import json
import numpy as np
import ROOT

def get_results(fit_name: str, output_path: str, significance: float, no_line: int):
    json_file = os.path.join(output_path, f"histo.{fit_name}.json")

    with open(json_file, 'r') as file:
        my_json = json.load(file)

    if f"intensity{no_line}_in_cts" in my_json[fit_name]["fit_parameters"]["line"].keys():
        gm =            my_json[fit_name]["fit_parameters"]["line"][f"intensity{no_line}_in_cts"]["mode"] 
        E0_counts_L68 = my_json[fit_name]["fit_parameters"]["line"][f"intensity{no_line}_in_cts"]["range_min"] 
        E0_counts_U68 = my_json[fit_name]["fit_parameters"]["line"][f"intensity{no_line}_in_cts"]["range_max"] 
        E0_counts_U90 = my_json[fit_name]["fit_parameters"]["line"][f"intensity{no_line}_in_cts"]["upper_limit"]

        root_file = f"histo_marginalized.{fit_name}.root"
        file_root = ROOT.TFile.Open(os.path.join(output_path, root_file), "READ")
        if not os.path.join(output_path, root_file):
            print(f"Error - ROOT file does not exist for {fit_name} - exit here.")
            return None, None, None

        h = file_root.Get(f"h1_histo_fitter_model_parameter_intensity{no_line}")
        h_clone = file_root.Get(f"h1_histo_fitter_model_parameter_intensity{no_line}")

        binmax = h.GetMaximumBin()
        for k in range(0,h.GetNbinsX(),1):
            h.SetBinContent(k, h.GetBinContent(k+h.FindBin(gm)-1)) # shift to left (global mode centred on 0)
        p = [0, 0.3415, 0.46, 0.50, 0.84, 0.90, 0.95]
        q = np.zeros(len(p), dtype=float)
        h_clone.GetQuantiles(len(p), q, np.array(p, dtype=float)) # GetQuantiles(n_division, quantiles, probsum[%])
        qt_34 = q[1] - q[0]

        # delete histos and close file
        h = None
        h_clone = None
        file_root.Close()

        # compare global mode with 34% quantile
        if gm/qt_34 < significance:
            return E0_counts_U90, 0, 0
        else:
            return gm, E0_counts_L68, E0_counts_U68

