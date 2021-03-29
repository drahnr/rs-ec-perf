/// Inlined functions for implementing GF arithmetics for SSSE3 instruction sets.
///

#ifndef _GF16_SSE_H_
#define _GF16_SSE_H_


include!("gf16_tables.h");

// SSE2 #include <emmintrin.h>

// SSSE3 #include <tmmintrin.h>






#[inline(always)]
fn linear_transform_8x8_128b( tab_l: __m128i, tab_h: __m128i, v: __m128i, mask_f: __m128i ) -> __m128i
{
	_mm_shuffle_epi8(tab_l,v&mask_f)^_mm_shuffle_epi8(tab_h,_mm_srli_epi16(v,4)&mask_f)
}



//////////////  GF(16)  /////////////////////////////




#[inline(always)]
fn tbl_gf16_squ( a: __m128i ) -> __m128i
{
	let tab_l: __m128i = _mm_load_si128((__m128i const *) __gf16_squ );
	_mm_shuffle_epi8(tab_l,a)
}


#[inline(always)]
fn tbl_gf16_squ_x8( a: __m128i ) -> __m128i
{
	let tab_l: __m128i = _mm_load_si128((__m128i const *) __gf16_squ_x8 );
	_mm_shuffle_epi8(tab_l,a)
}


#[inline(always)]
fn tbl_gf16_inv( a: __m128i ) -> __m128i
{
	let tab_l: __m128i = _mm_load_si128((__m128i const *) __gf16_inv );
	_mm_shuffle_epi8(tab_l,a);
}

#[inline(always)]
fn tbl_gf16_log( a: __m128i ) -> __m128i
{
	let tab_l: __m128i = _mm_load_si128((__m128i const *) __gf16_log );
	_mm_shuffle_epi8(tab_l,a)
}

#[inline(always)]
fn tbl_gf16_exp( a: __m128i ) -> __m128i
{
	let tab_l: __m128i = _mm_load_si128((__m128i const *) __gf16_exp );
	_mm_shuffle_epi8(tab_l,a)
}

#[inline(always)]
fn tbl_gf16_mul_const( unsigned char a , __m128i b ) -> __m128i
{
	let tab_l: __m128i = _mm_load_si128((__m128i const *) (__gf256_mul+  ((unsigned)a)*32 ));
	_mm_shuffle_epi8(tab_l,b)
}



#[inline(always)]
fn tbl_gf16_mul_0x8( __m128i b ) -> __m128i
{
	let tab_l: __m128i = _mm_load_si128((__m128i const *) (__gf256_mul+  8*32 ));
	_mm_shuffle_epi8(tab_l,b)
}




#[inline(always)]
fn tbl_gf16_mul( __m128i a , __m128i b ) -> __m128i
{
	let mask_f: __m128i = _mm_load_si128((__m128i const *) __mask_low);
	let log_16: __m128i = _mm_load_si128((__m128i const *) __gf16_log);
	let exp_16: __m128i = _mm_load_si128((__m128i const *) __gf16_exp);

	let la: __m128i = _mm_shuffle_epi8(log_16,a);
	let lb: __m128i = _mm_shuffle_epi8(log_16,b);
	let la_lb: __m128i = _mm_add_epi8(la,lb);

	let r0: __m128i = _mm_shuffle_epi8(exp_16, _mm_sub_epi8(la_lb, mask_f&_mm_cmpgt_epi8(la_lb,mask_f) ) );
	return r0;
}



#[inline(always)]
fn tbl_gf16_mul_log( __m128i a , __m128i logb , __m128i mask_f ) -> __m128i
{
	let la: __m128i = tbl_gf16_log( a );
	let la_lb: __m128i = _mm_add_epi8(la,logb);
	tbl_gf16_exp( _mm_sub_epi8(la_lb, mask_f&_mm_cmpgt_epi8(la_lb,mask_f) ) )
}

#[inline(always)]
fn tbl_gf16_mul_log_log( __m128i loga , __m128i logb , __m128i mask_f ) -> __m128i
{
	let la_lb: __m128i = _mm_add_epi8(loga,logb);
	tbl_gf16_exp( _mm_sub_epi8(la_lb, mask_f&_mm_cmpgt_epi8(la_lb,mask_f) ) )
}





/////////////////////////////  GF(256) ////////////////////////////////////////




#[inline(always)]
fn tbl_gf256_mul_const( a: u8 , b: __m128i ) -> __m128i
{
	let mask_f: __m128i = _mm_load_si128((__m128i const *) __mask_low);
	let tab_l: __m128i = _mm_load_si128((__m128i const *) (__gf256_mul+  ((unsigned)a)*32 ));
	let tab_h: __m128i = _mm_load_si128((__m128i const *) (__gf256_mul+  ((unsigned)a)*32 + 16 ));
	linear_transform_8x8_128b( tab_l , tab_h , b , mask_f )
}


/// use log table
#[inline(always)]
fn tbl_gf256_mul( a: __m128i , b: __m128i ) -> __m128i
{
	let mask_f: __m128i = _mm_load_si128((__m128i const *) __mask_low);
	let log_16: __m128i = _mm_load_si128((__m128i const *) __gf16_log);
	let exp_16: __m128i = _mm_load_si128((__m128i const *) __gf16_exp);

	let a0: __m128i = a&mask_f;
	let a1: __m128i = _mm_srli_epi16(a,4)&mask_f;
	let b0: __m128i = b&mask_f;
	let b1: __m128i = _mm_srli_epi16(b,4)&mask_f;

	let la0: __m128i = _mm_shuffle_epi8(log_16,a0);
	let la1: __m128i = _mm_shuffle_epi8(log_16,a1);
	let lb0: __m128i = _mm_shuffle_epi8(log_16,b0);
	let lb1: __m128i = _mm_shuffle_epi8(log_16,b1);

	let la0b0: __m128i = _mm_add_epi8(la0,lb0);
	let la1b0: __m128i = _mm_add_epi8(la1,lb0);
	let la0b1: __m128i = _mm_add_epi8(la0,lb1);
	let la1b1: __m128i = _mm_add_epi8(la1,lb1);

	let r0: __m128i = _mm_shuffle_epi8(exp_16, _mm_sub_epi8(la0b0, mask_f&_mm_cmpgt_epi8(la0b0,mask_f) ) );
	let r1: __m128i = _mm_shuffle_epi8(exp_16, _mm_sub_epi8(la1b0, mask_f&_mm_cmpgt_epi8(la1b0,mask_f) ) )
			^_mm_shuffle_epi8(exp_16, _mm_sub_epi8(la0b1, mask_f&_mm_cmpgt_epi8(la0b1,mask_f) ) );
	let r2: __m128i = _mm_shuffle_epi8(exp_16, _mm_sub_epi8(la1b1, mask_f&_mm_cmpgt_epi8(la1b1,mask_f) ) );

	_mm_slli_epi16(r1^r2,4)^r0^tbl_gf16_mul_0x8(r2)
}


/*
#[inline(always)]
fn tbl_gf16_squ_sl4( a: __m128i ) -> __m128i
{
	let tab_l: __m128i = _mm_load_si128((__m128i const *) __gf16_squ_sl4 );
	_mm_shuffle_epi8(tab_l,a)
}

#[inline(always)]
fn tbl_gf256_squ( a: __m128i ) -> __m128i
{
	let mask_f: __m128i = _mm_load_si128((__m128i const *) __mask_low);
	let a0: __m128i = a&mask_f;
	let a1: __m128i = _mm_srli_epi16(a,4)&mask_f;
	let a0squ: __m128i = tbl_gf16_squ(a0);
	let a1squ_sl4: __m128i = tbl_gf16_squ_sl4(a1);
	let a1squ_x8: __m128i = tbl_gf16_squ_x8( a1 );
	a1squ_sl4^a0squ^a1squ_x8
}
*/

#[inline(always)]
fn __m128i tbl_gf256_inv( a: __m128i )
{
    if true {
        // faster
    	let mask_f: __m128i = _mm_load_si128((__m128i const *) __mask_low);
    	let a0: __m128i = a&mask_f;
    	let a1: __m128i = _mm_srli_epi16(a,4)&mask_f;
    	let a0_a1: __m128i = a0^a1;
    	let a1squx8: __m128i = tbl_gf16_squ_x8( a1 );
    	let a0xa0_a1: __m128i = tbl_gf16_mul( a0 , a0_a1 );
    	let denominator: __m128i = a1squx8^a0xa0_a1;
    	let _denominator: __m128i = tbl_gf16_inv( denominator );
    	let b1: __m128i = tbl_gf16_mul( _denominator , a1 );
    	let a1inv: __m128i = tbl_gf16_inv(a1);
    	let b01: __m128i = tbl_gf16_mul( a0_a1 , a1inv );
    	let b01: __m128i = tbl_gf16_mul( b01 , b1 );
    	let a1x8: __m128i = tbl_gf16_mul_0x8( a1 );
    	let a0inv: __m128i = tbl_gf16_inv(a0);
    	let a1x8xb1_1: __m128i = tbl_gf16_mul( a1x8 , b1 ) ^ _mm_set1_epi8(1);
    	let b02: __m128i = tbl_gf16_mul( a0inv , a1x8xb1_1 );
    	let b0: __m128i = _mm_setzero_si128();
    	b0 |= _mm_andnot_si128( _mm_cmpeq_epi8(b0,a1inv), b01 ) |  _mm_andnot_si128( _mm_cmpeq_epi8(b0,a0inv), b02 );
    	_mm_slli_epi16(b1,4)^b0
    } else {
        // slow
    	let a2: __m128i = tbl_gf256_squ(a);
    	let a3: __m128i = tbl_gf256_mul(a2,a);
    	let a6: __m128i = tbl_gf256_squ(a3);
    	let a7: __m128i = tbl_gf256_mul(a6,a);
    	let ae: __m128i = tbl_gf256_squ(a7);
    	let af: __m128i = tbl_gf256_mul(ae,a);
    	let af1: __m128i = tbl_gf256_squ(af);
    	let af2: __m128i = tbl_gf256_squ(af1);
    	let af3: __m128i = tbl_gf256_squ(af2);
    	let a7f: __m128i = tbl_gf256_mul(af3,a7);
    	tbl_gf256_squ(a7f)
    }
}




#[inline(always)]
fn tbl_gf256_set_value( a: u8 ) -> __m128i { _mm_set1_epi8(a) }


#[inline(always)]
fn _tbl_gf256_set_value( unsigned char * b, a: u8 ) {
	_mm_storeu_si128( (__m128i *)b , _mm_set1_epi8(a) );
}


#[inline(always)]
fn tbl_gf256_get_1st_value( a: __m128i ) -> u8 { _mm_extract_epi16(a,0) & 0xff }


