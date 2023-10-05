# How to run the code
From terminal, outside any legend sif container, type

```bash
    $ ./run.sh config/config.json
```


# Config entries
Put here some useful paths and info to retrieve data:
```bash
    {
        "gamma-src-code": "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis", // where you have located the folder
        "output": "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/output/2023-09-04/", // where you want to store any fit output
        "prodenv": "/lfs/l1/legend/data/public/prodenv/prod-blind/ref/", // folder containing processed data
        "version": "v01.06", // version of processed data
```

Below, you need to provide input for which data you want to inspect (eg. periods/runs/detector type/cut)
```bash        
        "dataset": {
            "periods": ["p03", "p04"], // list of periods to inspect
            "runs": [["r000", "r001"], ["r000"]], // for each period, list of runs to inspect
            "status": ["on", "no_psd"], // detectors status to keep in the analysis
            "detectors": "BEGe", // choose among 'BEGe', 'COAX', 'PPC', 'ICPC', 'All', 'single'
            "cut": "raw", // choose among 'raw', 'LArC', 'LArAC'
```
In particular, you can build your specific dataset in the following way:
- _golden_ dataset: `"status": ["on"]`
- _silver_ dataset: `"status": ["on", "no_psd"]`
- _bronze_ dataset: `"status": ["on", "no_psd", "ac"]`
For each dataset, a different set of cuts is applied over flags saved in the skimmed files.

If you select `"detectors": "All"`, you'll inspect the spectrum given by the sum of counts coming from all detectors. 
If you select `"detectors": "singe"`, you'll inspect the spectrum of each detector separately. However, the analysis here considers only the K lines because of their higher statistics (useful to perform studies on the spatial distribution of K-40 and K-42 contaminants).

Finally, there is a specific entry for building ROOT files, ie histograms containing the spectra that you can inspect:
```bash
        
            "histogram": {
                "overwrite": true, 
                "folder": "raw_BEGe", // name of the folder where you want to store the ROOT files; it's created under 'src/root_files/'
                "bin-width": 1, // bin width in keV
                "x-min-keV": 0, // minimum in keV of the x-axis of the energy spectrum
                "x-max-keV": 5000 // maximum in keV of the x-axis of the energy spectrum
        }
    }
```
Select `"overwrite": true` if you want to re-generate the ROOT files (ie overwrite the previous ROOT folder, if exists), otherwise `"overwrite": false` if you want to skip the generation of ROOT files because already done in the past.

# `Src` content
Apart from the core of the code, here you can find:
- `gamma-fitter`: folder containing the BAT-C++ Bayesian analysis 
- `root_files`: folder containing the ROOT files that are produced starting from the processed skimmed files; if not present, this is created the first time you run the code
- `settings`: folder containing several tools and files produced along the running of the code (ex: exposure values, resolution values, list of available detectors)