# Changelog

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
