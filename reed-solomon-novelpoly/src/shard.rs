
use std::iter;
use crate::field::FieldAdd;
use crate::errors::*;

// pub trait Shard<F: FieldAdd>:
// Clone + AsRef<[u8]> + AsMut<[u8]> + AsMut<[[u8; <F as FieldAdd>::FIELD_BYTES]]> + AsRef<F::ElementAsBytes> + iter::FromIterator<F::ElementAsBytes> + From<Vec<u8>>
//     where [(); <F as FieldAdd>::FIELD_BYTES]: Sized

pub trait Shard<const FIELD_BYTES: usize>:
Clone + AsRef<[u8]> + AsMut<[u8]> + AsMut<[[u8; FIELD_BYTES]]> + AsRef<[[u8; FIELD_BYTES]]> + iter::FromIterator<[u8; FIELD_BYTES]> + From<Vec<u8>>
    where [u8; FIELD_BYTES]: Sized
{
  	type Inner;
  	fn into_inner(self) -> Self::Inner;

    fn set_chunk(&mut self, chunk_index: usize, chunk_data: [u8; FIELD_BYTES]) {
        AsMut::<[[u8; FIELD_BYTES]]>::as_mut(self)[chunk_index] = chunk_data;
    }

    fn get_chunk(&self, chunk_index: usize) -> [u8; FIELD_BYTES] {
        AsRef::<[[u8; FIELD_BYTES]]>::as_ref(&self)[chunk_index]
    }    
}

impl<T, const FIELD_BYTES: usize> Shard<FIELD_BYTES> for T
where
	T: Clone
		+ AsRef<[u8]>
		+ AsMut<[u8]>
		+ AsMut<[[u8; FIELD_BYTES]]>
		+ AsRef<[[u8; FIELD_BYTES]]>
		+ iter::FromIterator<[u8; FIELD_BYTES]>
	+ From<Vec<u8>>
    , [u8; FIELD_BYTES]: Sized
    
{
	type Inner = Self;
	fn into_inner(self) -> Self::Inner {
		self
	}

}




