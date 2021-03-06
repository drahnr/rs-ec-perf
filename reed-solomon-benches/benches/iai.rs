use iai::black_box;
use reed_solomon_benches::novelpoly::WrappedShard;
use reed_solomon_benches::*;
use reed_solomon_tester::*;

const N_SHARDS_MANY: usize = 2000;

fn bench_iai_novelpoly_roundtrip() {
	roundtrip(
		novelpoly::encode::<WrappedShard>,
		novelpoly::reconstruct::<WrappedShard>,
		black_box(BYTES),
		N_SHARDS_MANY,
	)
	.unwrap();
}

fn bench_iai_novelpoly_encode() {
	novelpoly::encode::<WrappedShard>(black_box(BYTES), N_SHARDS_MANY).unwrap();
}

iai::main!(bench_iai_novelpoly_roundtrip, bench_iai_novelpoly_encode);
