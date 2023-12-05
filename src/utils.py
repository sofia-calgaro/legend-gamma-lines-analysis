import os
import json
import scipy
import numpy as np
import ROOT
from matplotlib import pyplot as plt

def get_results(fit_name: str, output_path: str, significance: float, no_line: int):
    """ Get best-fit and errors at 68%. If not, an upper limit at 90% CI is derived."""
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


def get_survival_fraction(fit_name: str, raw_output_path: str, no_line: int, N_samples: int, truncated_sf: bool, save_path=None):
    """Extract survival fractions starting from marginalized posterior pdfs."""
    raw_filename = os.path.join(raw_output_path, f"histo_marginalized.{fit_name}.root")
    raw_file_root = ROOT.TFile.Open(raw_filename, "READ")
    if not raw_file_root:
        print(f"Error - raw ROOT file does not exist for {fit_name} - exit here.")
        return None, None

    larC_filename = raw_filename.replace("raw", "LAr_C")
    if not os.path.isfile(larC_filename):
        larC_filename = raw_filename.replace("raw", "LArC")
        if not os.path.isfile(larC_filename):
             print(f"Error - LAr C ROOT file does not exist for {fit_name} - exit here.")
             return None, None
    larC_file_root = ROOT.TFile.Open(larC_filename, "READ")

    larAC_filename = raw_filename.replace("raw", "LAr_AC")
    if not os.path.isfile(larAC_filename):
        larAC_filename = raw_filename.replace("raw", "LArAC")
        if not os.path.isfile(larAC_filename):
             print(f"Error - LAr AC ROOT file does not exist for {fit_name} - exit here.")
             return None, None
    larAC_file_root = ROOT.TFile.Open(larAC_filename, "READ")

    # get posterior pdfs for intensities
    raw_posterior = raw_file_root.Get(f"h1_histo_fitter_model_parameter_intensity{no_line}")
    larC_posterior = larC_file_root.Get(f"h1_histo_fitter_model_parameter_intensity{no_line}")
    larAC_posterior = larAC_file_root.Get(f"h1_histo_fitter_model_parameter_intensity{no_line}")

    # fill with edges and bin contents
    raw_bins = [raw_posterior.GetXaxis().GetBinLowEdge(bin) for bin in range(1, raw_posterior.GetNbinsX() + 2)]
    raw_counts = []
    for i in range(0, raw_posterior.GetNbinsX()):
        raw_counts.append( raw_posterior.GetBinContent(i) )

    larC_bins = [larC_posterior.GetXaxis().GetBinLowEdge(bin) for bin in range(1, larC_posterior.GetNbinsX() + 2)]
    larC_counts = []
    for i in range(0, larC_posterior.GetNbinsX()):
        larC_counts.append( larC_posterior.GetBinContent(i) )

    larAC_bins = [larAC_posterior.GetXaxis().GetBinLowEdge(bin) for bin in range(1, larAC_posterior.GetNbinsX() + 2)]
    larAC_counts = []
    for i in range(0, larAC_posterior.GetNbinsX()):
        larAC_counts.append( larAC_posterior.GetBinContent(i) )

    # generate distributions given the posterior histograms
    raw_sampled = scipy.stats.rv_histogram([raw_counts, raw_bins])
    larC_sampled = scipy.stats.rv_histogram([larC_counts, larC_bins])
    larAC_sampled = scipy.stats.rv_histogram([larAC_counts, larAC_bins])

    # sample N times 
    larC_survival_pdf = []
    larAC_survival_pdf = []
    for i in range(0, N_samples, 1):
        ratio_C = larC_sampled.rvs()/raw_sampled.rvs()
        ratio_AC = larAC_sampled.rvs()/raw_sampled.rvs()
        
        # truncated pdf
        if truncated_sf:
            if (ratio_C >= 0) & (ratio_C <= 1):
                larC_survival_pdf.append(ratio_C)
            if (ratio_AC >= 0) & (ratio_AC <= 1):
                larAC_survival_pdf.append(ratio_AC)
        # non-truncated pdf
        else:
            larC_survival_pdf.append(ratio_C)
            larAC_survival_pdf.append(ratio_AC)

    # define histos wit new pdfs
    surv_f_larC_hist  = plt.hist(larC_survival_pdf, bins = larC_posterior.GetNbinsX())
    surv_f_larAC_hist = plt.hist(larAC_survival_pdf, bins = larAC_posterior.GetNbinsX())

    # get mode
    surv_f_larC_mode = surv_f_larC_hist[1][np.where(surv_f_larC_hist[0] == np.max(surv_f_larC_hist[0]))]
    surv_f_larAC_mode = surv_f_larAC_hist[1][np.where(surv_f_larAC_hist[0] == np.max(surv_f_larAC_hist[0]))]

    # calculate central quantiles wrt global mode
    percentiles = [34.1, 84.1]
    
    larC_quantiles = np.percentile(larC_survival_pdf, percentiles)
    larAC_quantiles = np.percentile(larAC_survival_pdf, percentiles)
    """
    larC_distances_to_mode = np.abs(larC_survival_pdf - surv_f_larC_mode)
    larAC_distances_to_mode = np.abs(larAC_survival_pdf - surv_f_larAC_mode)
    larC_quantiles = [np.percentile(larC_distances_to_mode, 16) + surv_f_larC_mode, np.percentile(larC_distances_to_mode, 84) + surv_f_larC_mode]
    larAC_quantiles = [np.percentile(larAC_distances_to_mode, 16) + surv_f_larAC_mode, np.percentile(larAC_distances_to_mode, 84) + surv_f_larAC_mode]
    larC_quantiles = [k[0] for k in larC_quantiles]
    larAC_quantiles = [k[0] for k in larAC_quantiles]
    """
    """
    larC_low_lim, larC_upp_lim = interval(surv_f_larC_hist[1], surv_f_larC_hist[0], surv_f_larC_mode[0], 0.68)
    larC_quantiles = [larC_low_lim, larC_upp_lim]
    larAC_low_lim, larAC_upp_lim = interval(surv_f_larAC_hist[1], surv_f_larAC_hist[0], surv_f_larAC_mode[0], 0.68)
    larAC_quantiles = [larAC_low_lim, larAC_upp_lim]
    print(larC_quantiles)
    print(larAC_quantiles)
    print(surv_f_larC_mode)
    """

    # check if quantiles overcome 0 and 1 thresholds for non-truncated pdfs
    """
    if not truncated_sf:
        if larC_quantiles[0]<0:
            larC_quantiles[0] = 0
        if larC_quantiles[-1]>1:
            larC_quantiles[-1] = 1
        if larAC_quantiles[0]<0:
            larAC_quantiles[0] = 0
        if larAC_quantiles[-1]>1:
            larAC_quantiles[-1] = 1
    """

    # display quantiles
    for p, q in zip(percentiles, larC_quantiles):
        print(f"LArC --> {p}% percentile: {q}")
    for p, q in zip(percentiles, larAC_quantiles):
        print(f"LArAC --> {p}% percentile: {q}")

    # save posteriors in a given location (if specified)
    if save_path is not None:
        if not os.path.isdir(save_path) and not os.path.exists(save_path):
            os.makedirs(save_path)
        save_path = os.path.join(save_path, f"{fit_name}.json")
        my_dict = {
                "raw": {
                    "raw_counts": raw_counts,
                    "raw_bins": raw_bins,
                    },
                "LArC":{
                    "larC_counts": larC_counts,
                    "larC_bins": larC_bins,
                    "larC_survival_pdf": larC_survival_pdf,
                    "larC_quantiles": larC_quantiles
                    },
                "LArAC":{
                    "larAC_counts": larAC_counts,
                    "larAC_bins": larAC_bins,
                    "larAC_survival_pdf": larAC_survival_pdf,
                    "larAC_quantiles": larAC_quantiles
                    }
        }
        with open(save_path, 'w') as file:
            json.dump(my_dict, file, indent=4)
        print(f"Results saved in {save_path}")

        save_plot_path = save_path.replace("json", "pdf")  
        plt.figure(figsize=(10,6))
        plt.hist(larC_survival_pdf, bins = larC_posterior.GetNbinsX(), histtype = "step", facecolor="None", linewidth=1, color = "dodgerblue", label="LAr_C / raw")
        plt.hist(larAC_survival_pdf, bins = larAC_posterior.GetNbinsX(), histtype = "step", facecolor="None", linewidth=1, color = "limegreen", label="LAr_AC / raw")
        plt.axvline(x=larC_quantiles[0], color = "darkblue", linestyle="--", label="+-68%")
        plt.axvline(x=larC_quantiles[-1], color = "darkblue", linestyle="--")
        plt.axvline(x=larAC_quantiles[0], color = "darkblue", linestyle="--")
        plt.axvline(x=larAC_quantiles[-1], color = "darkblue", linestyle="--")
        plt.axvline(x=surv_f_larC_mode, color = "red", label="Global mode")
        plt.axvline(x=surv_f_larAC_mode, color = "red")
        plt.legend()
        plt.ylabel("Counts")
        plt.xlabel("Survival fraction")
        plt.tight_layout()
        plt.savefig(save_plot_path)
        print(f"Plot saved in {save_plot_path}")

    return [surv_f_larC_mode, larC_quantiles[0], larC_quantiles[-1]], [surv_f_larAC_mode, larAC_quantiles[0], larAC_quantiles[-1]]
