
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

    fn set_chunk(mut self, chunk_index: usize, chunk_data: [u8; FIELD_BYTES]) {
        AsMut::<[[u8; FIELD_BYTES]]>::as_mut(&mut self)[chunk_index] = chunk_data;
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

/// It deals with a slice of shards and verify that they are in sane 
/// condition to be used in ReedSolomon reconstruction
/// 
pub trait ShardHold<'a, S, const FIELD_BYTES: usize> :
Clone +
    AsRef<[Option<S>]> + AsMut<[Option<S>]>
where
    S: Shard<FIELD_BYTES>,

{
    ///Verifies if all shards have the same length and they can 
    ///be propely converted to a slice of underlying field elements
    ///return the uniform shard length
    fn verify_reconstructiblity(&self) -> Result<usize> {
        //if all shards empty there is nothig to reconstruct hence reject.
        let maybe_first_available_shard = AsRef::<[Option<S>]>::as_ref(&self).iter().find(|optional_shard| match optional_shard { Some(_) => true, None => false});
        let first_available_shard = match maybe_first_available_shard.as_ref() {
            None => Err(Error::PayloadSizeIsZero)?,
            Some(first_available_shard) => first_available_shard.as_ref().expect("Already has checked it is not none. q.e.d"),
        };
        
        let uniform_shard_len = AsRef::<[u8]>::as_ref(&first_available_shard).len();

        if uniform_shard_len == 0 {
            Err(Error::ZeroLengthShards)?;
        }

        if uniform_shard_len % FIELD_BYTES != 0 {
            Err(Error::UndivisableShardLength {shard_length: uniform_shard_len, field_bytes: FIELD_BYTES} )?;            
        }
            
        for optional_shard in self.as_ref() {//  { AsRef::<[Option<S>]>.as_ref(&self)
            match optional_shard {
                Some(shard) => if  AsRef::<[u8]>::as_ref(shard).len() !=  uniform_shard_len {
                    Err(Error::InconsistentShardLengths{ first: uniform_shard_len, other: AsRef::<[u8]>::as_ref(shard).len() })?;
                }
                _ => (),
            }
        }

        return Ok(uniform_shard_len);
    }

    ///make set of the shard to have exactly as many shard as
    ///the number of symbols in an encoded word, by either adding
    ///empty shards or removing extra shards.
    fn equalize_shards_number_with_code_block_length(self: &'a mut Self, code_block_length: usize) -> Vec<Option<&'a S>> {
        let gap = code_block_length.saturating_sub(self.as_ref().len()); //== max(code_block_length - self.as_ref().len(), 0): minimum number of missing shards, some received shard might be None

        //This might be too naive you might be removing none empty shards and leaving empty shards in place. nonetheless given the placement of the shard in the slice are important it is not possible to rescue beyond block length data without major restructuring of the reconstruction code
			AsMut::<[Option<S>]>::as_mut(self).iter().map(|s| s.as_ref()).take(code_block_length).chain(std::iter::repeat(None).take(gap)).collect::<Vec<_>>()
    }

}


