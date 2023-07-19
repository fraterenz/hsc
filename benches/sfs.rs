use criterion::black_box;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::BenchmarkId;
use criterion::Criterion;
use hsc::process::NeutralMutationPoisson;
use hsc::stemcell::{Sfs, StemCell};
use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use rand_distr::Poisson;

fn from_neutral_rate(c: &mut Criterion) {
    let mut group = c.benchmark_group("from_neutral_rate");
    let mutations = vec![1usize, 2usize, 3usize];
    let cells = [
        StemCell::with_set_of_mutations(mutations.clone()),
        StemCell::with_set_of_mutations(mutations),
    ];
    let verbosity = 0;

    for rate in [10., 1., 0.1, 0.01, 0.001, 0.0001].iter() {
        group.bench_with_input(BenchmarkId::from_parameter(rate), rate, |b, _| {
            b.iter(|| {
                let mut rng = ChaChaRng::seed_from_u64(26);
                let distributions = NeutralMutationPoisson(Poisson::new(*rate).unwrap());
                Sfs::from_cells(
                    black_box(&cells),
                    black_box(&distributions),
                    verbosity,
                    &mut rng,
                );
            })
        });
    }

    group.finish();
}

criterion_group!(benches, from_neutral_rate);
criterion_main!(benches);
