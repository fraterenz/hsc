## Output
The output of the simulations are saved in the directory given by the user.
If the directory path given by the user is `/path/to/save/`, then the output of the simulations will look like this:
```$ tree /path/to/save/
/path/to/save/
test/
├── 10cells
│   ├── burden # single cell mutational burden mapping
│   │   ├── 1 # the more recent timepoint
│   │   │   ├── 0.json  # run 0
│   │   │   ├── 1.json  # run 1
│   │   │   └── ...
│   │   ├── 2
│   │   │   ├── 0.json
│   │   │   ├── 1.json
│   │   │   └── ...
│   │   └── ...
│   ├── sfs # site frequency spectrum vector
│   │   ├── 1
│   │   │   ├── 0.json
│   │   │   ├── 1.json
│   │   │   └── ...
│   │   ├── 2
│   │   │   ├── 0.json
│   │   │   ├── 1.json
│   │   │   └── ...
│   │   └── ...
│   └── variant_fraction # the subclones' abbundance
│       ├── 1
│       │   ├── 0.csv
│       │   ├── 1.csv
│       │   └── ...
│       ├── 2
│       │   ├── 0.csv
│       │   ├── 1.csv
│       │   └── ...
│       └── ...
└── rates
    ├── 0.csv
    ├── 1.csv
    └── ...
```
that is those measurements are saved at the end of the simulations, for each timepoint, for each run.

Note that the order of the timepoints is reversed, hence timepoint 1 is the one saved at last.

### Measurements
- **burden**: a json file with keys being the number of mutations (single-cell mutational burden x-axis) and values being the cells with those mutations (single-cell mutational burden y-axis)
- **sfs:** a json file with keys being the jcells (x-axis) and values being the number of variants with jcells (y-axis).
- **variant_fraction**: the abbundance of all subclones
- **rates:** the birth-rates of the subclones, where the first entry represents the birth-rate of the wild-type `b0`

