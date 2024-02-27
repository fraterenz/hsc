# Changelog
The semantic versioning is kind of random.

## 3.0.4
### BugFix
- When cloning a cell upon division, **do not** set `last_division_t` of the daughter cell to 0

## 3.0.3
### BugFix
- Fix bug with times and background mutations
- When cloning a cell upon division, set `last_division_t` of the daughter cell to 0
- Fix when restarting, set timer to 0 and set all `last_division_t` to 0

## 3.0.2
### BugFix
- Fix bug for background mutations, when the interdivision time is small (<0.01) returns `None`

## 3.0.1
### BugFix
- forgot to commit new module `src/proliferation.rs`

## 3.0.0
- implement the background mutations for the exponential growing phase, by default it's on now.
- remove all the remaining symmetric stuff

### BugFix
- update all the background mutations while saving (even without taking a snapshot)

## 2.2.14
- increase the range of sigma up to 1

## 2.2.13
- in the first exponential phase, start with one cell with one neutral mutation

## 2.2.12
- increase fitness range up to 4

## 2.2.11
### BugFix
- when subsampling save the whole population as well with new flag `--save-population`
- fix some default parameters in cli

## 2.2.10
- increase range std to 0.5

## 2.2.9
- increase the max number of subclones from 1000 to 1200
- increase range fitness and std up to 2 and 0.2

## 2.2.8
- increase the max number of subclones from 800 to 1000

## 2.2.7
- increase the max number of subclones from 750 to 800
- decrease the min possible range for `std` to 0.001

## 2.2.6
- decrease the max number of subclones from 1000 to 750 such that jobs with smaller tau can hopefully finish in 1h.

## 2.2.5
- increase the max number of subclones from 400 to 1000.

- increase the max number of subclones from 400 to 1000.

## 2.2.4
- update dep with dependabot

## 2.2.3
- update dep with dependabot

## 2.2.2
- remove all the asymmetric stuff
### BugFix
- remove factor 2 in the neutral mutation rate

## 2.2.1
### BugFix
- fix subsampling after the exponential phase and also better subsampling when we subsample with different number of cells but at the same time. Note that subsampling after the exponential phase occurs after one extra division from the Moran phase.
- fix skipping the first timepoint to subsample

## 2.2.0
- allow subsampling at specific timepoints with different number of cells.

## 2.1.1
### BugFix
- at the end of the simulation, update the background mutations for all cells.

## 2.1.0
- increase the number of clones from 200 to 400

## 2.0.0
- change the behaviour of leaving empty arguments: empty parameters were sampled from intervals before (for ABC), but now instead they are constant and not sampled from any interval.

## 1.4.0
- remove parameter `r` and replace it by `--tau`, being `1/r`
### BugFix
- fix the rate for the mutational burden for `r` different than one

## 1.3.5
- rename the proliferation rate of the wild type from `mu0` to `r`.
- generate random `r` when not provided between range of 0.1 and 5.

## 1.3.4
- change random ranges of parameters mean: \[0.01, 0.4\], std: \[0.005, 0.1\] and mu: \[0.1, 20\].
- increase the number of clones to 200

## 1.3.3
- ensure the mean and the std are not too high (< 1 and smaller than 0.1 respectively)
- increase the max number of iterations

## 1.3.2
### BugFix
- fix bug in the variant fraction

## 1.3.1
- increase the number of clones to 160
- decrease the range of mu
### BugFix
- fix saving the variants with subsampling (issue #58)

## 1.3.0
- save timepoints using the `time` from the process instead of the lenght of the timepoints

## 1.2.0
- `clap_app` arguments `nb_snapshots` and `snapshots`
- `clippy`

## 1.1.2
- Increase the number of clones to 130
### BugFix
- fix filename for the std

## 1.1.1
- Generate random parameters for `u`, `mean` and `std` if not provided.

## 1.1.0
- Save files using the filename with the rate of fit mutants per cell division (divide by two since we are in symmetric case), mean and std of the fitness distribution.

## 1.0.0
- Remove the clever implementation and just simulate every cell and every mutations (with uuid).

## 0.19.0
- simulations start at year 0 and end at year `year` not `year+1`
### BugFix
-fix the two mutational rates: the two parameters of the Moran and Exponential growing phases were mixed up

## 0.18.0
- the first timepoint is saved at year 0 not at year 1

## 0.17.0
- clap app has two neutral rates, one for the Moran process and the other optional one for the exponential growing phase
### BugFix
- the neutral rate for the exponential growth was used also for the Moran process

## 0.16.1
- reduce the neutral mutation rate during the growing phase by setting it to `neutral_rate / 4`

## 0.16.0
- neutral scenario in clap app
- initial exponential growing phase

## 0.15.0
- sample the fitness coefficients from the gamma distribution

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
