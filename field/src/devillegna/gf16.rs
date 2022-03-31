/// Arithmetic in GF(16) and GF(256)

// gf4 := gf2[x]/x^2+x+1
#[inline(always)]
fn gf4_mul_2(a: u8) -> u8 {
    let r = a << 1;
    r ^= (a >> 1) * 7;
    r
}

#[inline(always)]
fn gf4_mul_3(a: u8) -> u8 {
    let msk = (a - 2) >> 1;
    (msk & (a * 3)) | ((~msk) & (a - 1))
}

#[inline(always)]
fn gf4_mul(a: u8, b: u8) -> u8 {
    let r = a * (b & 1);
    r ^ (gf4_mul_2(a) * (b >> 1))
}

#[inline(always)]
fn gf4_squ(a: u8) -> u8 {
    a ^ (a >> 1)
}

#[inline(always)]
fn gf4_inv(a: u8) -> u8 {
    a ^ (a >> 1)
}

///////

#[inline(always)]
fn gf4v_mul_2_u32(a: u32) -> u32 {
    let bit0 = a & 0x55555555;
    let bit1 = a & 0xaaaaaaaa;
    (bit0 << 1) ^ bit1 ^ (bit1 >> 1)
}

#[inline(always)]
fn gf4v_mul_3_u32(a: u32) -> u32 {
    let bit0 = a & 0x55555555;
    let bit1 = a & 0xaaaaaaaa;
    (bit0 << 1) ^ bit0 ^ (bit1 >> 1)
}

#[inline(always)]
fn gf4v_mul_u32(a: u32, b: u8) -> u32 {
    let bit0_b = ((u32) 0) - ((u32)(b & 1));
    let bit1_b = ((u32) 0) - ((u32)((b >> 1) & 1));
    (a & bit0_b) ^ (bit1_b & gf4v_mul_2_u32(a))
}

#[inline(always)]
fn _gf4v_mul_u32_u32(u32 a0, u32 a1, u32 b0, u32 b1) -> u32 {
    let c0 = a0 & b0;
    let c2 = a1 & b1;
    let c1_ = (a0 ^ a1) & (b0 ^ b1);
    ((c1_ ^ c0) << 1) ^ c0 ^ c2
}

#[inline(always)]
fn gf4v_mul_u32_u32(a: u32, b: u32) -> u32 {
    let a0 = a & 0x55555555;
    let a1 = (a >> 1) & 0x55555555;
    let b0 = b & 0x55555555;
    let b1 = (b >> 1) & 0x55555555;
    _gf4v_mul_u32_u32(a0, a1, b0, b1)
}

#[inline(always)]
fn gf4v_squ_u32(a: u32) -> u32 {
    let bit1 = a & 0xaaaaaaaa;
    a ^ (bit1 >> 1)
}

//////////////////////////////////////////////////////////////////////////////////

#[inline(always)]
fn gf16_is_nonzero(a: u8) -> u8 {
    unsigned a4 = a & 0xf;
    unsigned r = ((unsigned) 0) - a4;
    r >>= 4;
    r & 1
}

// gf16 := gf4[y]/y^2+y+x
#[inline(always)]
fn gf16_mul(a: u8, b: u8) -> u8 {
    let a0 = a & 3;
    let a1 = (a >> 2);
    let b0 = b & 3;
    let b1 = (b >> 2);
    let a0b0 = gf4_mul(a0, b0);
    let a1b1 = gf4_mul(a1, b1);
    let a0b1_a1b0 = gf4_mul(a0 ^ a1, b0 ^ b1) ^ a0b0 ^ a1b1;
    let a1b1_x2 = gf4_mul_2(a1b1);
    ((a0b1_a1b0 ^ a1b1) << 2) ^ a0b0 ^ a1b1_x2
}

#[inline(always)]
fn gf16_squ(a: u8) -> u8 {
    let a0 = a & 3;
    let a1 = (a >> 2);
    let a1 = gf4_squ(a1);
    let a1squ_x2 = gf4_mul_2(a1);
    (a1 << 2) ^ a1squ_x2 ^ gf4_squ(a0)
}

#[inline(always)]
fn gf16_inv(a: u8) -> u8 {
    let a2 = gf16_squ(a);
    let a4 = gf16_squ(a2);
    let a8 = gf16_squ(a4);
    let a6 = gf16_mul(a4, a2);
    gf16_mul(a8, a6)
}

#[inline(always)]
fn gf16_mul_4(a: u8) -> u8 {
    (((a << 2) ^ a) & (8 + 4)) ^ gf4_mul_2(a >> 2)
}

#[inline(always)]
fn gf16_mul_8(a: u8) -> u8 {
    let a0 = a & 3;
    let a1 = a >> 2;
    (gf4_mul_2(a0 ^ a1) << 2) | gf4_mul_3(a1)
}

////////////

// gf16 := gf4[y]/y^2+y+x

#[inline(always)]
fn gf16v_mul_u32(a: u32, b: u8) -> u32 {
    let axb0 = gf4v_mul_u32(a, b);
    let axb1 = gf4v_mul_u32(a, b >> 2);
    let a0b1 = (axb1 << 2) & 0xcccccccc;
    let a1b1 = axb1 & 0xcccccccc;
    let a1b1_2 = a1b1 >> 2;
    axb0 ^ a0b1 ^ a1b1 ^ gf4v_mul_2_u32(a1b1_2)
}

#[inline(always)]
fn _gf16v_mul_u32_u32(u32 a0, u32 a1, u32 a2, u32 a3, u32 b0, u32 b1, u32 b2, u32 b3) -> u32 {
    let c0 = _gf4v_mul_u32_u32(a0, a1, b0, b1);
    let c1_ = _gf4v_mul_u32_u32(a0 ^ a2, a1 ^ a3, b0 ^ b2, b1 ^ b3);

    let c2_0 = a2 & b2;
    let c2_2 = a3 & b3;
    let c2_1_ = (a2 ^ a3) & (b2 ^ b3);
    let c2_r0 = c2_0 ^ c2_2;
    let c2_r1 = c2_0 ^ c2_1_;
    //u32 c2 = c2_r0^(c2_r1<<1);
    // GF(4) x2: (bit0<<1)^bit1^(bit1>>1);
    ((c1_ ^ c0) << 2) ^ c0 ^ (c2_r0 << 1) ^ c2_r1 ^ (c2_r1 << 1)
}

#[inline(always)]
fn gf16v_mul_u32_u32(a: u32, b: u32) -> u32 {
    let a0 = a & 0x11111111;
    let a1 = (a >> 1) & 0x11111111;
    let a2 = (a >> 2) & 0x11111111;
    let a3 = (a >> 3) & 0x11111111;
    let b0 = b & 0x11111111;
    let b1 = (b >> 1) & 0x11111111;
    let b2 = (b >> 2) & 0x11111111;
    let b3 = (b >> 3) & 0x11111111;

    _gf16v_mul_u32_u32(a0, a1, a2, a3, b0, b1, b2, b3)
}

#[inline(always)]
fn gf256v_reduce_u32(a: u32) -> u8 {
    // https://godbolt.org/z/7hirMb
    let aa = a.to_le_bytes();
    aa[0] ^ aa[1] ^ aa[2] ^ aa[3]
}

#[inline(always)]
fn gf16v_reduce_u32(a: u32) -> u8 {
    let r256 = gf256v_reduce_u32(a);
    (r256 & 0xf) ^ (r256 >> 4)
}

#[inline(always)]
fn gf16v_squ_u32(a: u32) -> u32 {
    let a2 = gf4v_squ_u32(a);
    a2 ^ gf4v_mul_2_u32((a2 >> 2) & 0x33333333)
}

#[inline(always)]
fn gf16v_mul_4_u32(a: u32) -> u32 {
    let a1 = a & 0xcccccccc;
    let a0 = (a << 2) & 0xcccccccc;
    a0 ^ a1 ^ gf4v_mul_2_u32(a1 >> 2)
}

#[inline(always)]
fn gf16v_mul_8_u32(a: u32) -> u32 {
    let a1 = a & 0xcccccccc;
    let a0 = (a << 2) & 0xcccccccc;
    gf4v_mul_2_u32(a0 ^ a1) | gf4v_mul_3_u32(a1 >> 2)
}

////////////

#[inline(always)]
fn gf256_is_nonzero(a: u8) -> u8 {
    unsigned a8 = a;
    unsigned r = ((unsigned) 0) - a8;
    r >>= 8;
    r & 1
}

// gf256 := gf16[X]/X^2+X+xy
#[inline(always)]
fn gf256_mul(a: u8, b: u8) -> u8 {
    let a0 = a & 15;
    let a1 = (a >> 4);
    let b0 = b & 15;
    let b1 = (b >> 4);
    let a0b0 = gf16_mul(a0, b0);
    let a1b1 = gf16_mul(a1, b1);
    let a0b1_a1b0 = gf16_mul(a0 ^ a1, b0 ^ b1) ^ a0b0 ^ a1b1;
    let a1b1_x8 = gf16_mul_8(a1b1);
    ((a0b1_a1b0 ^ a1b1) << 4) ^ a0b0 ^ a1b1_x8
}

#[inline(always)]
fn gf256_mul_gf16(a: u8, gf16_b: u8) -> u8 {
    let a0 = a & 15;
    let a1 = (a >> 4);
    let b0 = gf16_b & 15;
    let a0b0 = gf16_mul(a0, b0);
    let a1b0 = gf16_mul(a1, b0);
    a0b0 ^ (a1b0 << 4)
}

#[inline(always)]
fn gf256_squ(a: u8) -> u8 {
    let a0 = a & 15;
    let a1 = (a >> 4);
    let a1 = gf16_squ(a1);
    let a1squ_x8 = gf16_mul_8(a1);
    (a1 << 4) ^ a1squ_x8 ^ gf16_squ(a0)
}

#[inline(always)]
fn gf256_inv(a: u8) -> u8 {
    // 128+64+32+16+8+4+2 = 254
    let a2 = gf256_squ(a);
    let a4 = gf256_squ(a2);
    let a8 = gf256_squ(a4);
    let a4_2 = gf256_mul(a4, a2);
    let a8_4_2 = gf256_mul(a4_2, a8);
    let a64_ = gf256_squ(a8_4_2);
    let a64_ = gf256_squ(a64_);
    let a64_ = gf256_squ(a64_);
    let a64_2 = gf256_mul(a64_, a8_4_2);
    let a128_ = gf256_squ(a64_2);
    gf256_mul(a2, a128_)
}

////////

#[inline(always)]
fn gf256v_mul_u32(a: u32, b: u8) -> u32 {
    let axb0 = gf16v_mul_u32(a, b);
    let axb1 = gf16v_mul_u32(a, b >> 4);
    let a0b1 = (axb1 << 4) & 0xf0f0f0f0;
    let a1b1 = axb1 & 0xf0f0f0f0;
    let a1b1_4 = a1b1 >> 4;

    axb0 ^ a0b1 ^ a1b1 ^ gf16v_mul_8_u32(a1b1_4)
}

#[inline(always)]
fn gf256v_squ_u32(a: u32) -> u32 {
    let a2 = gf16v_squ_u32(a);
    let ar = (a2 >> 4) & 0x0f0f0f0f;

    a2 ^ gf16v_mul_8_u32(ar)
}

#[inline(always)]
fn gf256v_mul_gf16_u32(a: u32, gf16_b: u8) -> u32 {
    gf16v_mul_u32(a, gf16_b)
}


// gf256 := gf16[X]/X^2+X+xy
#[inline(always)]
fn gf256v_mul_0x10_u32(a: u32) -> u32 {
    let a0 = a&0x0f0f0f0f;
    let a1 = a&0xf0f0f0f0;
    let a1x = gf16v_mul_8_u32(a1>>4);
    (a0<<4)^a1^a1x
}

#[inline(always)]
fn gf256v_mul_0x20_u32(a: u32) -> u32 {
    let a0 = gf4v_mul_2_u32(a&0x0f0f0f0f);
    let a1 = gf4v_mul_2_u32(a&0xf0f0f0f0);
    let a1x = gf16v_mul_8_u32(a1>>4);
    (a0<<4)^a1^a1x
}

#[inline(always)]
fn gf256v_mul_0x40_u32(a: u32) -> u32 {
    let a0 = gf16v_mul_4_u32(a&0x0f0f0f0f);
    let a1 = gf16v_mul_4_u32(a&0xf0f0f0f0);
    let a1x = gf16v_mul_8_u32(a1>>4);
    (a0<<4)^a1^a1x
}

#[inline(always)]
fn gf256v_mul_0x80_u32(a: u32) -> u32 {
    let a0 = gf16v_mul_8_u32(a&0x0f0f0f0f);
    let a1 = gf16v_mul_8_u32(a&0xf0f0f0f0);
    let a1x = gf16v_mul_8_u32(a1>>4);
    (a0<<4)^a1^a1x
}

