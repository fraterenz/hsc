## Output
The output of the simulations are saved in the directory given by the user.
If the directory path given by the user is `/path/to/save/`, then the output of the simulations will look like this:
```tree /path/to/save/
/path/to/save/
├── sfs
│   ├── 1 # <- timepoint 1, the oldest
│   │   └── 0.json # <- run 0
│   │   └── 1.json # <- run 1
│   │   └── ...
│   ├── 2
│   │   └── 0.json
│   │   └── 1.json
│   │   └── ...
│   ├── ...
├── sfs_neutral
│   ├── 1
│   │   └── 0.json
│   │   └── 1.json
│   │   └── ...
│   ├── 2
│   │   └── 0.json
│   │   └── 1.json
│   │   └── ...
│   ├── ...
└── variant_fraction
    ├── 1
    │   └── 0.csv
    │   └── 1.csv
    │   └── ...
    ├── 2
    │   └── 0.csv
    │   └── 1.csv
    │   └── ...
    ├── ...
```
that is 3 measurements are saved at the end of the simulations, for each timepoint, for each run.
Note that the order of the timepoints is reversed, hence timepoint 1 is the one saved at last.

### Measurements
- **sfs**: a json file with keys being the jcells (sfs x-axis) and values being the mutations present in jcells (sfs y-axis)
- **sfs_neutral**: a json file with keys being the jcells (sfs x-axis) and values being the mutations present in jcells (sfs y-axis) for the wild-type cells only, i.e. the cells being in the subclone with id 0, that is with no fitness advantage
- **variant_fraction**: the proportion of cells in all subclones

