use std::iter;
use crate::field::FieldAdd;

//pub trait Shard<F: FieldAdd>:
// Clone + AsRef<[u8]> + AsMut<[u8]> + AsMut<[[u8; <F as FieldAdd>::FIELD_BYTES]]> + AsRef<F::Element> + iter::FromIterator<F::Element> + From<Vec<u8>>
// //    where [u8; <F as FieldAdd>::FIELD_BYTES]:  [(); Sized]
//     where [(); <F as FieldAdd>::FIELD_BYTES]: Sized

pub trait Shard<F: FieldAdd>:
Clone + AsRef<[u8]> + AsMut<[u8]> + AsMut<[F::Element]> + AsRef<F::Element> + iter::FromIterator<F::Element> + From<Vec<u8>>
//    where F::Element : Sized
{
  	type Inner;
  	fn into_inner(self) -> Self::Inner;

}

// impl<T, F> Shard<F> for T
// where
// 	T: Clone
// 		+ AsRef<[u8]>
// 		+ AsMut<[u8]>
// 		+ AsMut<[[u8; 2]]>
// 		+ AsRef<[[u8; 2]]>
// 		+ iter::FromIterator<[u8; 2]>
// 	+ From<Vec<u8>>,
//     F: FieldAdd
// {
// 	type Inner = Self;
// 	fn into_inner(self) -> Self::Inner {
// 		self
// 	}
// }
