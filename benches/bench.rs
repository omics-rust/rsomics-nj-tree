use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_nj_tree(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-nj-tree");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let tsv = manifest.join("tests/golden/dist.tsv");
    c.bench_function("rsomics-nj-tree golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .arg(tsv.to_str().unwrap())
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_nj_tree);
criterion_main!(benches);
