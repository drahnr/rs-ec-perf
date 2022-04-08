// A shard with a even number of elements, which can sliced into 2 byte haps
use crate::FieldAdd;
use std::marker::PhantomData;

#[derive(Clone, Debug)]
pub struct WrappedShard<F: FieldAdd> {
    inner: Vec<u8>,
    _marker: PhantomData<*const F>
    
}

impl<F: FieldAdd> WrappedShard<F> {
    /// Wrap `data`.
    pub fn new(mut data: Vec<u8>) -> Self {
        if data.len() & 0x01 == 0x01 {
            data.push(0);
        }

        WrappedShard::<F> { inner: data, _marker: PhantomData }
    }

    /// Unwrap and yield inner data.
    pub fn into_inner(self) -> Vec<u8> {
        self.inner
    }
}

impl<F: FieldAdd> From<Vec<u8>> for WrappedShard<F> {
    fn from(data: Vec<u8>) -> Self {
        Self::new(data)
    }
}

impl<F: FieldAdd> AsRef<[u8]> for WrappedShard<F> {
    fn as_ref(&self) -> &[u8] {
        self.inner.as_ref()
    }
}

impl<F: FieldAdd> AsMut<[u8]> for WrappedShard<F> {
    fn as_mut(&mut self) -> &mut [u8] {
        self.inner.as_mut()
    }
}

impl<F: FieldAdd> AsRef<[[u8; F::FIELD_BYTES]]> for WrappedShard<F> {
    fn as_ref(&self) -> &[[u8; F::FIELD_BYTES]] {
        assert_eq!(self.inner.len() & 0x01, 0);
        if self.inner.is_empty() {
            return &[];
        }
        unsafe { ::std::slice::from_raw_parts(&self.inner[0] as *const _ as _, self.inner.len() / 2) }
    }
}

impl<F: FieldAdd> AsMut<[[u8; F::FIELD_BYTES]]> for WrappedShard<F> {
    fn as_mut(&mut self) -> &mut [[u8; F::FIELD_BYTES]] {
        let len = self.inner.len();
        assert_eq!(len & 0x01, 0);

        if self.inner.is_empty() {
            return &mut [];
        }
        unsafe { ::std::slice::from_raw_parts_mut(&mut self.inner[0] as *mut _ as _, len / 2) }
    }
}


impl<F: FieldAdd> std::iter::FromIterator<[u8; F::FIELD_BYTES]>  for WrappedShard<F>
{

    fn from_iter<T>(iterable: T) -> Self        
        where  T: IntoIterator<Item = [u8; F::FIELD_BYTES]>
    {
         let iter = iterable.into_iter();
         let (l, _) = iter.size_hint();
        let mut inner = Vec::with_capacity(F::FIELD_BYTES);

        for cur_chunk in iter {
            cur_chunk.iter().map(|a| inner.push(*a));
        }

        debug_assert_eq!(inner.len() & 0x01, 0);
        WrappedShard::<F> { inner: inner, _marker: PhantomData }
    }
}
