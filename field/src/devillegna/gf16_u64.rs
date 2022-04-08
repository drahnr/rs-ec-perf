/// Inlined functions for implementing GF arithmetics for u64 arch.

include!("gf16.h")



#[inline(always)]
fn gf4v_mul_2_u64( u64 a ) -> u64
{
	u64 bit0 = a&0x5555555555555555ull;
	u64 bit1 = a&0xaaaaaaaaaaaaaaaaull;
	(bit0<<1)^bit1^(bit1>>1)
}


#[inline(always)]
fn gf4v_mul_3_u64( u64 a ) -> u64
{
	u64 bit0 = a&0x5555555555555555ull;
	u64 bit1 = a&0xaaaaaaaaaaaaaaaaull;
	(bit0<<1)^bit0^(bit1>>1)
}


#[inline(always)]
fn gf4v_mul_u64( u64 a , b: u8 ) -> u64
{
	u64 bit0_b = ((u64)0)-((u64)(b&1));
	u64 bit1_b = ((u64)0)-((u64)((b>>1)&1));
	(a&bit0_b)^(bit1_b&gf4v_mul_2_u64(a))
}


#[inline(always)]
fn _gf4v_mul_u64_u64( u64 a0 , u64 a1 , u64 b0 , u64 b1 ) -> u64
{
	u64 c0 = a0&b0;
	u64 c2 = a1&b1;
	u64 c1_ = (a0^a1)&(b0^b1);
	((c1_^c0)<<1)^c0^c2
}

#[inline(always)]
fn gf4v_mul_u64_u64( u64 a , u64 b ) -> u64
{
	u64 a0 = a&0xaaaaaaaaaaaaaaaaull;
	u64 a1 = (a>>1)&0xaaaaaaaaaaaaaaaaull;
	u64 b0 = b&0xaaaaaaaaaaaaaaaaull;
	u64 b1 = (b>>1)&0xaaaaaaaaaaaaaaaaull;

	_gf4v_mul_u64_u64( a0 , a1 , b0 , b1 )
}



#[inline(always)]
fn gf4v_squ_u64( u64 a ) -> u64
{
	u64 bit1 = a&0xaaaaaaaaaaaaaaaaull;
	a^(bit1>>1)
}



//////////////////////////////////////////////////////////////////////////////////


// gf16 := gf4[y]/y^2+y+x


#[inline(always)]
fn gf16v_mul_u64( u64 a , unsigned char b ) -> u64
{
	u64 axb0 = gf4v_mul_u64( a , b );
	u64 axb1 = gf4v_mul_u64( a , b>>2 );
	u64 a0b1 = (axb1<<2)&0xccccccccccccccccull;
	u64 a1b1 = axb1&0xccccccccccccccccull;
	u64 a1b1_2 = a1b1 >>2;

	axb0 ^ a0b1 ^ a1b1 ^ gf4v_mul_2_u64( a1b1_2 )
}



#[inline(always)]
fn _gf16v_mul_u64_u64( u64 a0 , u64 a1 , u64 a2 , u64 a3 , u64 b0 , u64 b1 , u64 b2 , u64 b3 ) -> u64
{
	u64 c0 = _gf4v_mul_u64_u64( a0 , a1 , b0 , b1 );
	u64 c1_ =  _gf4v_mul_u64_u64( a0^a2 , a1^a3 , b0^b2 , b1^b3 );

	u64 c2_0 = a2&b2;
	u64 c2_2 = a3&b3;
	u64 c2_1_ = (a2^a3)&(b2^b3);
	u64 c2_r0 = c2_0^c2_2;
	u64 c2_r1 = c2_0^c2_1_;
	//u64 c2 = c2_r0^(c2_r1<<1);
	// GF(4) x2: (bit0<<1)^bit1^(bit1>>1);
	((c1_^c0)<<2) ^c0^ (c2_r0<<1)^c2_r1^(c2_r1<<1)
}

#[inline(always)]
fn gf16v_mul_u64_u64( u64 a , u64 b ) -> u64
{
	u64 a0 = a&0x1111111111111111ull;
	u64 a1 = (a>>1)&0x1111111111111111ull;
	u64 a2 = (a>>2)&0x1111111111111111ull;
	u64 a3 = (a>>3)&0x1111111111111111ull;
	u64 b0 = b&0x1111111111111111ull;
	u64 b1 = (b>>1)&0x1111111111111111ull;
	u64 b2 = (b>>2)&0x1111111111111111ull;
	u64 b3 = (b>>3)&0x1111111111111111ull;

	_gf16v_mul_u64_u64( a0 , a1 , a2 , a3 , b0 , b1 , b2 , b3 )
}

#[inline(always)]
fn gf256v_reduce_u64( u64 a ) -> u8
{
	uint32_t * aa = (uint32_t *)(&a);
	uint32_t r = aa[0]^aa[1];
	gf256v_reduce_u32( r )
}

#[inline(always)]
fn gf16v_reduce_u64( u64 a ) -> u8
{
	uint8_t r256 = gf256v_reduce_u64( a );
	(r256&0xf)^(r256>>4);
}



#[inline(always)]
fn gf16v_squ_u64( u64 a ) -> u64
{
	let a2 = gf4v_squ_u64( a );
	a2 ^ gf4v_mul_2_u64( (a2>>2)& 0x3333333333333333ull )
}

#[inline(always)]
fn gf16v_mul_8_u64( u64 a ) -> u64
{
	let a1 = a&0xccccccccccccccccull;
	let a0 = (a<<2)&0xccccccccccccccccull;
	gf4v_mul_2_u64(a0^a1)|gf4v_mul_3_u64(a1>>2)
}





//////////////////////////////////////////////////////////


#[inline(always)]
fn gf256v_mul_u64( u64 a , b: u8 ) -> u64
{
	u64 axb0 = gf16v_mul_u64( a , b );
	u64 axb1 = gf16v_mul_u64( a , b>>4 );
	u64 a0b1 = (axb1<<4) & 0xf0f0f0f0f0f0f0f0ull;
	u64 a1b1 = axb1&0xf0f0f0f0f0f0f0f0ull;
	u64 a1b1_4 = a1b1 >>4;

	axb0 ^ a0b1 ^ a1b1 ^ gf16v_mul_8_u64( a1b1_4 )
}


#[inline(always)]
fn gf256v_squ_u64( u64 a ) -> u64
{
	u64 a2 = gf16v_squ_u64( a );
	u64 ar = (a2>>4)&0x0f0f0f0f0f0f0f0full;

	a2 ^ gf16v_mul_8_u64( ar )
}


#[inline(always)]
fn gf256v_mul_gf16_u64( a: u64, gf16_b: u8 ) -> u64
{
	gf16v_mul_u64( a , gf16_b )
}



