use std::iter;
use std::fmt::Debug;
use crate::field::FieldAdd;
use crate::errors::*;

pub trait Shard<F: FieldAdd>:
    Clone
    + AsRef<[u8]>
    + AsMut<[u8]>
    + AsMut<[[u8; F::FIELD_BYTES]]>
    + AsRef<[[u8; F::FIELD_BYTES]]>
    + iter::FromIterator<[u8; F::FIELD_BYTES]>
    + From<Vec<u8>>
    + Debug
where
    [u8; F::FIELD_BYTES]: Sized,
{
    type Inner;
    fn into_inner(self) -> Self::Inner;

    fn set_chunk(&mut self, chunk_index: usize, chunk_data: [u8; F::FIELD_BYTES]) {
        AsMut::<[[u8; F::FIELD_BYTES]]>::as_mut(self)[chunk_index] = chunk_data;
    }

    fn get_chunk(&self, chunk_index: usize) -> [u8; F::FIELD_BYTES] {
        AsRef::<[[u8; F::FIELD_BYTES]]>::as_ref(&self)[chunk_index]
    }
}

impl<T, F> Shard<F> for T
where
    T: Clone
        + AsRef<[u8]>
        + AsMut<[u8]>
        + AsMut<[[u8; F::FIELD_BYTES]]>
        + AsRef<[[u8; F::FIELD_BYTES]]>
        + iter::FromIterator<[u8; F::FIELD_BYTES]>
    + From<Vec<u8>>
    + Debug
    ,
    F: FieldAdd,
    [u8; F::FIELD_BYTES]: Sized,
{
    type Inner = Self;
    fn into_inner(self) -> Self::Inner {
        self
    }
}
