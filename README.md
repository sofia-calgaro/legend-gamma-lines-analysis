# How to run the code
From terminal, outside any legend sif container, type

```bash
    $ ./run.sh config/config.json
```


# Config entries
Put here some useful paths:
```bash
    {
        "gamma-src-code": "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis", // where you have located the folder
        "output": "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/output/2023-09-04/", // where you want to store any fit output
        "histo-folder": "/lfs/l1/legend/users/calgaro/legend_analysis/legend-gamma-lines-analysis/src/root_files/0_25keV-golden/", // where you can find the generated ROOT files
```

Below, you need to provide input for data to read and data to inspect (eg. periods/runs/detector type/cut)
```bash        
        "dataset": {
            "prodenv": "/lfs/l1/legend/data/public/prodenv/prod-blind/ref/",
            "version": "v01.06",
            "comments": ["The user can set the list of periods, runs and detector types ('single','All','BEGe', 'ICPC', 'COAX','PPC' )",
                        "If they want to study only some detector then set a list [] in the 'detectors' key",
                        "The user can choose the cut to apply among ['raw',  'LAr AC', 'LAr C']"], 
            "periods": ["p03", "p04"],
            "runs": [["r000", "r001"], ["r000"]], // one list of runs for each period
            "detectors": "BEGe",
            "cut": "raw"
        },
```

Below, specify entries needed for exposure:
```bash
        "exposure": {
            "comments": ["The user can set the time unit of exposure evaluation ('sec', 'min', 'hour', 'day', 'yr')",
                    "The user can choose the status of detectors to include in the exposure evaluation among 'on', 'off', 'ac', 'no_psd' or a combination of them."],
            "time-unit": "yr", 
            "status": ["on"] // status of detector to take into account
        }
    }  
```
