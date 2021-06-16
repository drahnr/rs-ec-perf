
use std::iter;
use crate::field::FieldAdd;
use crate::errors::*;

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
    fn get_chunk(&self, chunk_index: usize) -> &[u8; <F as FieldAdd>::FIELD_BYTES];
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
        AsMut::<[F::ElementAsBytes]>::as_ref(&mut self)[chunk_index] = chunk_data;
    }

    fn get_chunk(&self, chunk_index: usize) -> &[u8; <F as FieldAdd>::FIELD_BYTES] {
        AsRef::<[F::ElementAsBytes]>::as_ref(&mut self)[chunk_index]
    }

    

}

/// It deals with a slice of shards and verify that they are in sane 
/// condition to be used in ReedSolomon reconstruction
/// 
pub trait ShardHold<S,F> :
    Clone + AsRef<[Option<S>]> where
    S: Shard<F>,
    F: FieldAdd,
    [u8; F::FIELD_BYTES]: Sized,


{
    ///Verifies if all shards have the same length and they can 
    ///be propely converted to a slice of underlying field elements
    ///return the uniform shard length
    fn verify_reconstructiblity(&self) -> Result<usize>;

    ///make set of the shard to have exactly as many shard as
    ///the number of symbols in an encoded word, by either adding
    ///empty shards or removing extra shards.
    fn equalize_shards_number_with_code_block_length(&self, code_block_length: usize) -> Self;
}

impl <S,F,T> ShardHold<S,F>  for T where
    T:
    Clone + AsRef<[Option<S>]>,
    S: Shard<F>,
    F: FieldAdd,
	[u8; F::FIELD_BYTES]: Sized,
    
{
    fn verify_reconstructiblity(&self) -> Result<usize>  {
        //if all shards empty there is nothig to reconstruct hence reject.
        let first_available_shard = self.into_iter().find_map(|optional_shard| optional_shard).unwrap_or(Err("No shard is available"));
        let uniform_shard_len = first_available_shard?.len();

        if uniform_shard_len == 0 {
            Err("Cannot reconstruct from length zero shards")?;
        }

        if uniform_shard_len % Shard::F::FIELD_BYTES == 0 {
            Err("Cannot properly divide shards into field elements")?;
        }
            
        for shard in self.as_ref() {//  { AsRef::<[Option<S>]>.as_ref(&self)
            match shard {
                Some(s) => match s.len() {
                    uniform_shard_len => (),
                    _ => Err(Error::InconsistentShardLengths{ first: uniform_shard_len, other: s.len })?,
                },
                _ => (),
            }
        }

        return uniform_shard_len;
    }

    fn equalize_shards_number_with_code_block_length(&self, code_block_length: usize) -> Self {
        let gap = code_block_length.saturating_sub(self.as_ref().len()); //minimum number of missing shards, some received shard might be None

        //This might be too naive you might be removing none empty shards and leaving empty shards in place. nonetheless given the placement of the shard in the slice are important it is not possible to rescue beyond block length data without major restructuring of the reconstruction code
		let rearrenged_shards =
			self.as_ref().into_iter().take(code_block_length).chain(std::iter::repeat(None).take(gap)).collect::<Vec<_>>();
    }

}

