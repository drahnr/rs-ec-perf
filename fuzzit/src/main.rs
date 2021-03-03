use honggfuzz::fuzz;

use rs::*;

fn main() {
	// Here you can parse `std::env::args and
	// setup / initialize your project

	novel_poly_basis::setup();

	// You have full control over the loop but
	// you're supposed to call `fuzz` ad vitam aeternam
	loop {
		// The fuzz macro gives an arbitrary object (see `arbitrary crate`)
		// to a closure-like block of code.
		// For performance reasons, it is recommended that you use the native type
		// `&[u8]` when possible.
		// Here, this slice will contain a "random" quantity of "random" data.
		fuzz!(|data: [u8; novel_poly_basis::K * 2]| {
			roundtrip(novel_poly_basis::encode, novel_poly_basis::reconstruct, &data);
		});
	}
}
