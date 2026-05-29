/// Utilities to save the Moran process.
///
/// We want to save at different timepoints (snapshots) some statistics or
/// quantities of a Moran process, the only process that can be saved for now.
use std::{
    fs,
    path::{Path, PathBuf},
};

use anyhow::Context;
use derive_builder::Builder;
use enumset::{EnumSet, EnumSetType};
use log::debug;
use sosa::ReactionRates;

use crate::{
    MAX_SUBCLONES,
    genotype::{MutationalBurden, Sfs, SingleCellMutations},
    stemcell::StemCell,
    subclone::{CloneId, SubClones, save_phylogeny, save_variant_fraction},
    write2file,
};

/// The statistics/measurements we want to save from the simulations.
#[derive(Hash, Debug, EnumSetType)]
pub enum Stats2Save {
    Burden,
    Sfs,
    VariantFraction,
    SingleCellMutations,
    VariantsPhylogeny,
    /// Per-clone Gillespie rates owned by the Moran process. Useful when
    /// rates evolve over time (e.g. `--multihits` resampling).
    Rates,
}

/// Specify whether to save the whole population or a subsampling.
#[derive(Debug, Clone, Copy)]
pub enum SavingCells {
    WholePopulation,
    Subsampling {
        sample_size: usize,
        population_as_well: bool,
    },
}

/// Specify when to save the measurements for a certain number of simulated cells.
#[derive(Debug, Clone)]
pub struct Snapshot {
    /// The number of cells to subsample
    pub cells2sample: usize,
    /// The time at which we subsample
    pub time: f32,
}

#[derive(Debug, Clone)]
pub struct SavingOptions {
    pub filename: PathBuf,
    pub save_population: bool,
}

#[derive(Debug, Copy, Clone, Builder)]
pub struct StatsConfig {
    #[builder(setter(custom), default)]
    enabled: EnumSet<Stats2Save>,
}

impl Default for StatsConfig {
    fn default() -> Self {
        let mut enabled = EnumSet::empty();
        enabled.insert(Stats2Save::Sfs);
        StatsConfig { enabled }
    }
}

impl StatsConfigBuilder {
    pub fn stats(&mut self, stats: impl IntoIterator<Item = Stats2Save>) -> &mut Self {
        let mut set = EnumSet::empty();
        for s in stats {
            set.insert(s);
        }
        self.enabled = Some(set);
        self
    }
}

fn make_path(
    path2dir: &Path,
    filename: &Path,
    tosave: Stats2Save,
    cells: usize,
    time: f32,
) -> anyhow::Result<PathBuf> {
    let mut timepoint = format!("{time:.1}").replace('.', "dot");
    timepoint.push_str("years");

    let path2dir = path2dir.join(format!("{cells}cells"));
    let path2file = match tosave {
        // only one file for the single cell mutational burden
        Stats2Save::SingleCellMutations => path2dir.join("mutations"),
        Stats2Save::VariantFraction => path2dir.join("variant_fraction"),
        Stats2Save::Burden => path2dir.join("burden"),
        Stats2Save::Sfs => path2dir.join("sfs"),
        Stats2Save::VariantsPhylogeny => path2dir.join("variant_phylogeny"),
        Stats2Save::Rates => path2dir.join("rates"),
    };
    let path2file = path2file.join(timepoint);

    fs::create_dir_all(&path2file).with_context(|| "Cannot create dir")?;
    debug!("creating dirs {path2file:#?}");

    Ok(path2file.join(filename))
}

pub(crate) fn save_it(
    path2dir: &Path,
    filename: &Path,
    time: f32,
    cells_with_idx: Vec<(&StemCell, CloneId)>,
    subclones: &SubClones,
    rates: &ReactionRates<MAX_SUBCLONES>,
    stats: &StatsConfig,
) -> anyhow::Result<()> {
    debug!("saving data at time {time}");
    let cells: Vec<&StemCell> = cells_with_idx.iter().map(|ele| ele.0).collect();
    let nb_cells = cells.len();

    debug!("saving {nb_cells} cells");

    if stats.enabled.contains(Stats2Save::Sfs) {
        if stats.enabled.contains(Stats2Save::SingleCellMutations) {
            let sc_mutations = SingleCellMutations::from_cells(&cells).unwrap_or_else(|_| {
                panic!("cannot create single cell mutations for the timepoint at time {time}")
            });
            sc_mutations.save(
                &make_path(
                    path2dir,
                    filename,
                    Stats2Save::SingleCellMutations,
                    nb_cells,
                    time,
                )?,
                time,
            )?;
            Sfs::from_sc_mutations(&sc_mutations)
                .unwrap_or_else(|_| panic!("cannot create SFS for timepoint at time {time}"))
                .save(
                    &make_path(path2dir, filename, Stats2Save::Sfs, nb_cells, time)?,
                    time,
                )?;
        } else {
            Sfs::from_cells(&cells)
                .unwrap_or_else(|_| panic!("cannot create SFS for timepoint at time {time}"))
                .save(
                    &make_path(path2dir, filename, Stats2Save::Sfs, nb_cells, time)?,
                    time,
                )?;
        }
    }

    if stats.enabled.contains(Stats2Save::SingleCellMutations) {
        SingleCellMutations::from_cells(&cells)
            .unwrap_or_else(|_| {
                panic!("cannot create single cell mutations for the timepoint at time {time}")
            })
            .save(
                &make_path(
                    path2dir,
                    filename,
                    Stats2Save::SingleCellMutations,
                    nb_cells,
                    time,
                )?,
                time,
            )?;
    }

    if stats.enabled.contains(Stats2Save::Burden) {
        MutationalBurden::from_cells(&cells)
            .unwrap_or_else(|_| panic!("cannot create burden for the timepoint at time {time}"))
            .save(
                &make_path(path2dir, filename, Stats2Save::Burden, nb_cells, time)?,
                time,
            )?;
    }

    if stats.enabled.contains(Stats2Save::VariantFraction) {
        save_variant_fraction(
            &SubClones::from(
                cells_with_idx
                    .into_iter()
                    .map(|(cell, id)| (cell.to_owned(), id))
                    .collect::<Vec<(StemCell, CloneId)>>(),
            ),
            &make_path(
                path2dir,
                filename,
                Stats2Save::VariantFraction,
                nb_cells,
                time,
            )?,
        )?;
    }

    if stats.enabled.contains(Stats2Save::VariantsPhylogeny) {
        save_phylogeny(
            subclones,
            &make_path(
                path2dir,
                filename,
                Stats2Save::VariantsPhylogeny,
                nb_cells,
                time,
            )?,
        )?;
    }

    if stats.enabled.contains(Stats2Save::Rates) {
        let path =
            make_path(path2dir, filename, Stats2Save::Rates, nb_cells, time)?.with_extension("csv");
        write2file(&rates.0, &path, None, false)
            .with_context(|| "cannot save the Gillespie rates of the subclones")?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MAX_SUBCLONES;

    fn unique_tmpdir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "hsc-test-{name}-{}-{:?}",
            std::process::id(),
            std::thread::current().id(),
        ));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn save_it_emits_phylogeny_when_enabled() {
        let dir = unique_tmpdir("save_it_phylogeny");
        let filename = PathBuf::from("0");
        let mut subclones = SubClones::new(vec![StemCell::new()], 4);
        subclones.get_mut_clone_unchecked(3).set_parent_id(1);

        let cells_with_idx = subclones.get_cells_with_clones_idx();
        let nb_cells = cells_with_idx.len();
        let rates = ReactionRates([0.0_f32; MAX_SUBCLONES]);

        let stats = StatsConfig {
            enabled: EnumSet::only(Stats2Save::VariantsPhylogeny),
        };

        save_it(
            &dir,
            &filename,
            0.0,
            cells_with_idx,
            &subclones,
            &rates,
            &stats,
        )
        .unwrap();

        let csv = dir
            .join(format!("{nb_cells}cells"))
            .join("variant_phylogeny")
            .join("0dot0years")
            .join("0.csv");
        assert!(csv.exists(), "expected {csv:?}");
        let content = std::fs::read_to_string(&csv).unwrap();
        let fields: Vec<&str> = content
            .strip_suffix(',')
            .unwrap_or(&content)
            .split(',')
            .collect();
        assert_eq!(fields.len(), MAX_SUBCLONES);
        assert_eq!(fields[0], "");
        assert_eq!(fields[3], "1");

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn save_it_emits_rates_when_enabled() {
        let dir = unique_tmpdir("save_it_rates");
        let filename = PathBuf::from("0");
        let subclones = SubClones::new(vec![StemCell::new(); 3], 4);
        let cells_with_idx = subclones.get_cells_with_clones_idx();
        let nb_cells = cells_with_idx.len();

        // Distinct entries at slot 0 (wild-type b0) and slot 1 (fit) so we can
        // assert the dense row layout end-to-end.
        let mut rates_arr = [0.5_f32; MAX_SUBCLONES];
        rates_arr[0] = 1.5;
        let rates = ReactionRates(rates_arr);

        let stats = StatsConfig {
            enabled: EnumSet::only(Stats2Save::Rates),
        };

        save_it(
            &dir,
            &filename,
            1.0,
            cells_with_idx,
            &subclones,
            &rates,
            &stats,
        )
        .unwrap();

        let csv = dir
            .join(format!("{nb_cells}cells"))
            .join("rates")
            .join("1dot0years")
            .join("0.csv");
        assert!(csv.exists(), "expected {csv:?}");
        let content = std::fs::read_to_string(&csv).unwrap();
        let fields: Vec<&str> = content
            .strip_suffix(',')
            .unwrap_or(&content)
            .split(',')
            .collect();
        assert_eq!(fields.len(), MAX_SUBCLONES);
        assert!((fields[0].parse::<f32>().unwrap() - 1.5).abs() < 1e-5);
        assert!((fields[1].parse::<f32>().unwrap() - 0.5).abs() < 1e-5);

        let _ = std::fs::remove_dir_all(&dir);
    }
}
