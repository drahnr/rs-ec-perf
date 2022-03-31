Arithmetic in GF(16) and GF(256) extracted from Ming-Shing Chen's work
https://github.com/fast-crypto-lab/rainbow-submission-round2/tree/master/Alternative_Implementation/ssse3/IIIc_Classic

But no SIMD code for GF(2^16) here.


The gf16.h has the code for multiplication in GF(16) and GF(256). The
fields are constructed as Eq2.11 of Ming-Shing Chen's thesis.

The gf16.c contains various static tables for GF(16). For example,
log/exp tables and multiplication tables.

The gf16_sse.h shows a SIMD implementation for both vector-vector and
vector-scalar multiplication for GF(16) and GF(256).
The vector-vector multiplication is implemented with log/exp tables of GF(16).


