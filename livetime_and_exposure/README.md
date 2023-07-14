Example scripts.

# Livetime
Store in a JSON file the livetime (in seconds) for each individual run and period.

# "get-exposure.py"
The macro returns and save a dictionary containing exposures (in kg*_unit_, where _unit_ can be specified by an user) channel by channel, run by run, period by period.

You can call the function from command line. To see available options use;
``` bash
$ python get_exposure --help
```

For instance,
``` bash
$ python get_exposure.py --livetime livetime_in_s.json --time_unit day --data  '{"p03": [0, 1, 2, 3, 4, 5],"p04": [0, 1, 2, 3, 6],"p05": [1, 2, 4],"p06": [0, 1, 2, 3]}' --status on
```
will get livetime values from "livetime_in_s.json", we express the exposure values in kg*day, we inspect p03/p04/p05/p06 runs, and we include only on detectors.

Or you can use the module from another script as
``` python
import get_exposure
exposures = get_exposure.main("livetime_in_s.json", "day", '{"p03": [0, 1, 2, 3, 4, 5],"p04": [0, 1, 2, 3, 6],"p05": [1, 2, 4],"p06": [0, 1, 2, 3]}', "on")
```
