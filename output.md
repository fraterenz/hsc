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
│   ├── genotypes # nested vec of cells with their proliferative events
│   │   ├── 1
│   │   │   ├── 0.json
│   │   │   ├── 1.json
│   │   │   └── ...
│   │   ├── 2
│   │   │   ├── 0.json
│   │   │   ├── 1.json
│   │   │   └── ...
│   │   └── ...
│   ├── sfs # site frequency spectrum vector
│   │   ├── 1
│   │   │   ├── 0.csv
│   │   │   ├── 1.csv
│   │   │   └── ...
│   │   ├── 2
│   │   │   ├── 0.csv
│   │   │   ├── 1.csv
│   │   │   └── ...
│   │   └── ...
│   ├── sfs_entropy # site frequency spectrum vector for the entropy
│   │   ├── 1
│   │   │   ├── 0.csv
│   │   │   ├── 1.csv
│   │   │   └── ...
│   │   ├── 2
│   │   │   ├── 0.csv
│   │   │   ├── 1.csv
│   │   │   └── ...
│   │   └── ...
│   ├── stats # site frequency spectrum vector for the entropy
│   │   │   ├── 0.csv
│   │   │   ├── 1.csv
│   │   │   └── ...
│   │   ├── 2
│   │   │   ├── 0.csv
│   │   │   ├── 1.csv
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
There is also a dir `stats` storing all the number of neutral mutations at the end of the simulation.
Note that the order of the timepoints is reversed, hence timepoint 1 is the one saved at last.

### Measurements
- **burden**: a json file with keys being the number of mutations (single-cell mutational burden x-axis) and values being the cells with those mutations (single-cell mutational burden y-axis)
- **genotype:** each entry represents a cell with its proliferative events (cell divisions)
- **sfs:** a vec where each entry represents a variant and the value stored indicates the number of cells with carrying that variant.
- **sfs_entropy:** same as sfs but the SFS is computed considering all cells at the timepoint of interest but counting only the variants that were present at a certain time in the past
- **stats:** a serialised struct storing the mapping between the proliferative events and the number of neutral mutations in the total population at the end of the simulation. Note that the entry `cell_count` is not correct, use this struct only with `poisson_mut_number`
- **variant_fraction**: the abbundance of all subclones
- **rates:** the birth-rates of the subclones, for now b0 * (1 + s) for all clones except wild-type which has proliferation rate equal to b0

