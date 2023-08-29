# Changelog
## 0.14.1
- simulations start at year zero and end at year `cli.years + 1`
- the first timepoint is saved at year 1 not at year 0

## 0.14.0
- rename clap app arg `cells2subsample` to `subsample` and allow for more than one entry for this arg
- move stats folder to the timepoint
- save the birth-rates of the subclones

## 0.13.1
- remove dbg print

## 0.13.0
- subsamples and change the dir for output

## 0.12.3
- SFS entropy by considering all cells at the timepoint of interest but counting only the variants that were present at a certain time in the past

## 0.12.2
- Replace the SFS for the proliferation events with the one with the variants
### BugFix
- Fix the single-cell mutational burden for the last timepoint by computing it `from_cells` and not `from_stats`

## 0.12.1
- Add the SFS for the proliferative events (not the variants).

## 0.12.0
### BugFix
- The SFS is actually the single cell mutational burden, so rename and refactor everything accordingly

## 0.11.2
### BugFix
- Fix the computation of the SFS

## 0.11.1
- do not raise an error when the genotype for the run cannot be found
- type changes: `nb_variants` is now `u64` and not `u8`; neutral mutations upon cell division are now `u16` not `u8`
- save `StatisticsMutations` into json

### BugFix
- fix the mutation rates by **not** dividing by 2 in the symmetric case. With v0.11.0 the proliferation scheme has changed, now in the symmetric case only one cell divides
- fix the neutral mutation rate by **not** dividing by the number of cells

## 0.11.0
- remove the neutral sfs

### BugFix
- fix the neutral rate
- fix the computation of the SFS for timepoints
- fix the acquisition of neutral mutations: we first cloned the proliferating cell and then mutated both cells. Now we just mutate one of the two cells, not both of them.

## 0.10.0
Start simulations at year 0 instead of 1.

### BugFix
- correct the Poisson rate for the neutral mutations.

## 0.9.0
Increase the number of clones from 40 to 60.

### BugFix
- now upon occurence of a fit mutation, assign cells to empty clones only (without overwriting older clones).

## 0.8.1
### Added
- Option `snapshots` to specify how many linespaced timepoints to save.

`sosa` version: 3.0.1

### BugFix
- Now save all timepoints

## 0.8.0
### Added
- When saving the timepoints, save also the sfs and the sfs neutral

## 0.7.3
### Added
- Benchmark sfs computation
- Swap the default hasher for sets and maps from the `std` with crate `rustc-hash`

### BugFix
Cell can have no neutral mutations upon cell-division: remove that nonzeroU.

## 0.7.2
Remove dependency `uuid` and `thiserror`, and change backend of the sfs computation.

## 0.7.1
Forgot to commit new files...

## 0.7.0
Rename files and fix bug in the SFS, still not sure whether it's correct or not.

## 0.6.0
### Added
- Simplify and disantangle the addition of neutral mutations and the positive selection clonal assignment.
- More doc
- New struct `NeutralMutationPoisson`

### BugFix
Fix bug in the neutral sfs (Fix #19).

## 0.5.0
Snapshots taken based on the time (not iterations). Many bugs fixed.

## 0.4.0
Tweek cli arguments.

## 0.3.0
Cli with `clap`.

## 0.2.0
Upload code.
