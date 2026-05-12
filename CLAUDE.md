# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

`hsc` is a Rust CLI tool for stochastic simulation of hematopoietic stem cell (HSC) dynamics during ageing. It models clonal competition between HSCs, tracking neutral mutations and positively selected (fit) mutations using the Gillespie algorithm via the [`sosa`](https://crates.io/crates/sosa) crate.

## Commands

```bash
# Build
cargo build --release

# Run all tests (unit + doc + quickcheck property tests)
cargo test --locked --all-features --all-targets
cargo test --locked --all-features --doc

# Run a single test
cargo test <test_name>

# Lint (CI enforces no warnings)
cargo clippy --all-targets --all-features -- -D warnings

# Format check
cargo fmt --all -- --check

# Format (apply)
cargo fmt

# Run the binary (example)
cargo run -- <DIR> moran
cargo run -- <DIR> exp-moran
RUST_LOG=info cargo run -- <DIR> moran   # verbose output
```

Log levels: `info`, `debug`, `trace` via `RUST_LOG`.

## Architecture

The simulation runs two types of stochastic processes, both driven by the Gillespie algorithm from `sosa`:

- **`Exponential`** (`process.rs`): exponential growth phase (embryonic development or post-treatment regrowth). Population size grows unbounded until a target cell count is reached.
- **`Moran`** (`process.rs`): fixed-size population phase (adult haematopoiesis). Every division is paired with a death to keep population constant. Snapshots are saved here.

Both implement `AdvanceStep<MAX_SUBCLONES>` from `sosa`, advancing one Gillespie step at a time.

### Typical simulation flow

1. Parse CLI args (`clap_app.rs`, `Cli::build()`) → `AppOptions`
2. Optionally run `Exponential` (dev/regrowth) then call `switch_to_moran()` to transition
3. Run `Moran` to the target age; save snapshots at configured timepoints
4. Optionally simulate treatment: subsample population (`into_subsampled`), regrow via `Exponential`, switch back to `Moran`
5. Write fitness rates to `<DIR>/rates/<idx>.csv`

### Key types

| Type | File | Role |
|---|---|---|
| `StemCell` | `stemcell.rs` | Individual cell; holds neutral `Mutation` UUIDs and last-division time |
| `SubClone` / `SubClones` | `subclone.rs` | Groups cells by fitness class; drives Gillespie rates |
| `Fitness` | `subclone.rs` | Enum: `Neutral`, `Fixed { s }`, `GammaSampled { shape, scale }` |
| `Distributions` | `subclone.rs` | Pre-computed Poisson/Bernoulli distributions for mutation sampling |
| `Proliferation` | `proliferation.rs` | Executes cell division; assigns neutral and fit mutations |
| `NeutralMutations` | `proliferation.rs` | `UponDivision` or `UponDivisionAndBackground` |
| `Snapshot` / `StatsConfig` | `snapshots.rs` | When and what to save (SFS, burden, variant fraction, mutations) |
| `Probs` / `ProbsPerYear` | `lib.rs` | Mutation rate parameters; converts per-year rates to per-division probabilities |

### Output format

See `output.md` for the full directory layout and file formats.

### Parallelism

Runs are parallelised with `rayon` by default. Use `--sequential` for sequential runs or `-d` (`--debug`) for a single debug run (10 cells, 20 iterations, high mutation rate).

Each run gets its own `ChaCha8Rng` stream derived from `--seed` + run index, making results reproducible.

### Constants

- `MAX_SUBCLONES = 3000` (`lib.rs`): hard upper bound on fit clones; increase if the simulation exits with an error about exceeding clone count.
- `TIME_AT_BIRTH = 9/12` years: delay between end of exponential phase and birth; used when assigning background mutations before the Moran phase.
