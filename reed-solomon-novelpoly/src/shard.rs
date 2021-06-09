use std::iter;
use crate::field::FieldAdd;

//pub trait Shard<F: FieldAdd>:
// Clone + AsRef<[u8]> + AsMut<[u8]> + AsMut<[[u8; <F as FieldAdd>::FIELD_BYTES]]> + AsRef<F::ElementAsBytes> + iter::FromIterator<F::ElementAsBytes> + From<Vec<u8>>
// //    where [u8; <F as FieldAdd>::FIELD_BYTES]:  [(); Sized]
//     where [(); <F as FieldAdd>::FIELD_BYTES]: Sized

pub trait Shard<F: FieldAdd>:
Clone + AsRef<[u8]> + AsMut<[u8]> + AsMut<[[u8; F::FIELD_BYTES]]> + AsRef<[[u8; F::FIELD_BYTES]]> + iter::FromIterator<[u8; F::FIELD_BYTES]> + From<Vec<u8>>
where
	[u8; F::FIELD_BYTES]: Sized,
{
  	type Inner;
  	fn into_inner(self) -> Self::Inner;

    fn set_chunk(&self, chunk_index: usize, chunk_data: &[u8; <F as FieldAdd>::FIELD_BYTES]);

}

impl<T, F> Shard<F> for T
where
	T: Clone
		+ AsRef<[u8]>
		+ AsMut<[u8]>
		+ AsMut<[[u8; F::FIELD_BYTES]]>
		+ AsRef<[[u8; F::FIELD_BYTES]]>
		+ iter::FromIterator<[u8; F::FIELD_BYTES]>
	+ From<Vec<u8>>,
    F: FieldAdd,
     [u8; F::FIELD_BYTES]: Sized
{
	type Inner = Self;
	fn into_inner(self) -> Self::Inner {
		self
	}

    fn set_chunk(&self, chunk_index: usize, chunk_data: &[u8; <F as FieldAdd>::FIELD_BYTES]) {
        AsMut::<[F::ElementAsBytes]>::as_ref(&mut self)[chunk_index] = chunk_data
    }
}
