use iai::black_box;
use reed_solomon_performance::*;

fn bench_iai_novel_poly_basis_roundtrip() {
	roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, black_box(BYTES), 2000).unwrap();
}

fn bench_iai_novel_poly_basis_encode() {
	novel_poly_basis::encode(black_box(BYTES), 2000).unwrap();
}

iai::main!(bench_iai_novel_poly_basis_roundtrip, bench_iai_novel_poly_basis_encode);
