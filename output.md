## Output
The output of the simulations are saved in the directory given by the user.
If the directory path given by the user is `/path/to/save/`, then the output of the simulations will look like this:
```$ tree /path/to/save/
/path/to/save/
test/
├── 10cells
│   ├── burden # single cell mutational burden mapping
│   │   ├── 0dot0years
│   │   │   ├── 0.json  # run 0
│   │   │   ├── 1.json  # run 1
│   │   │   └── ...
│   │   ├── 1dot0years
│   │   │   ├── 0.json
│   │   │   ├── 1.json
│   │   │   └── ...
│   │   └── ...
│   ├── mutations # single cell mut burden
│   │   ├── 0dot0years
│   │   │   ├── 0.parquet
│   │   │   ├── 1.parquet
│   │   │   └── ...
│   │   ├── 1dot0years
│   │   │   ├── 0.parquet
│   │   │   ├── 1.parquet
│   │   │   └── ...
│   │   └── ...
│   ├── sfs # site frequency spectrum vector
│   │   ├── 0dot0years
│   │   │   ├── 0.json
│   │   │   ├── 1.json
│   │   │   └── ...
│   │   ├── 1dot0years
│   │   │   ├── 0.json
│   │   │   ├── 1.json
│   │   │   └── ...
│   │   └── ...
│   ├── variant_fraction # the subclones' abbundance
│   │    ├── 0dot0years
│   │    │   ├── 0.csv
│   │    │   ├── 1.csv
│   │    │   └── ...
│   │    ├── 1dot0years
│   │    │   ├── 0.csv
│   │    │   ├── 1.csv
│   │    │   └── ...
│   │    └── ...
│   ├── variant_phylogeny # the subclones' parent_id
│   │    ├── 0dot0years
│   │    │   ├── 0.csv
│   │    │   ├── 1.csv
│   │    │   └── ...
│   │    ├── 1dot0years
│   │    │   ├── 0.csv
│   │    │   ├── 1.csv
│   │    │   └── ...
│   └    └── ...
└── rates
    ├── 0.csv
    ├── 1.csv
    └── ...
```
that is those measurements are saved at the end of the simulations, for each timepoint, for each run.

The file name of the runs (here `0.csv` and `1.csv` or `0.json` and `1.json`) have changed such that they incorportate the parameters that have been used to generate them, starting from version `v1.1.0`.

### Measurements
More info about the folders generated as the output of `hsc`:
- **`burden`**: this folder contains json files with representing the single-cell mutational burden, where keys are number of mutations found in cells (single-cell mutational burden x-axis) and values are the number of cells with these mutations (single-cell mutational burden y-axis),
- **`mutations`:** folder containing a parquet files representing the number of cells carrying the mutations present in the population,
- **`sfs`:** folder containing json files with keys being the jcells (x-axis) and values being the number of variants with jcells (y-axis),
- **`variant_fraction`**: folder containing csv files with the abbundance of all subclones, where the first entry represents the frequency of the wild-type clone with `b0` birth-rate,
- **`variant_phylogeny`**: folder containing csv files with the `parent_id` of all subclones — entry *i* is the id of the subclone that the cells in slot *i* descend from (set when a fit-mutation event placed them there).
The wild-type clone (slot 0) has `parent_id = None` because it has no parent. Every other slot defaults to `0` (descended from the wild-type clone) until a fit-mutation event reassigns it.
Note: as `parent_id` is a clone-level (not cell-level) property, when the snapshot is a subsample, `variant_phylogeny` always reflects the full-population lineage, and thus we can have `variant_phylogeny[i] != None` even when `variant_fractions[i] = 0`.
- **`rates`:** folder containing the csv files with the birth-rates of the subclones, where the first entry represents the birth-rate of the wild-type `b0`.

