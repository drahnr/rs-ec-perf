use iai::black_box;
use reed_solomon_performance::*;

fn bench_iai_novelpoly_roundtrip() {
	reed_solomon_tester::roundtrip(novelpoly::encode, novelpoly::reconstruct, black_box(BYTES), 2000).unwrap();
}

fn bench_iai_novelpoly_encode() {
	novelpoly::encode(black_box(BYTES), 2000).unwrap();
}

iai::main!(bench_iai_novelpoly_roundtrip, bench_iai_novelpoly_encode);
