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

use crate::{
    genotype::{MutationalBurden, Sfs, SingleCellMutations},
    stemcell::StemCell,
    subclone::{save_variant_fraction, SubClones},
};

/// The statistics/measurements we want to save from the simulations.
#[derive(Hash, Debug, EnumSetType)]
pub enum Stats2Save {
    Burden,
    Sfs,
    VariantFraction,
    SingleCellMutations,
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
    cells_with_idx: Vec<(&StemCell, usize)>,
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
                    .collect::<Vec<(StemCell, usize)>>(),
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

    Ok(())
}
