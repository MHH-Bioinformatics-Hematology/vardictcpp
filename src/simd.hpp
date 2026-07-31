#pragma once
// Portable SIMD helpers. Every kernel here has three backends selected at compile time from the
// architecture BASELINE instruction set, so no -march/-mavx flag is needed and the binary stays
// portable across any CPU of its architecture:
//   * x86-64  -> SSE2 (guaranteed on every x86-64 CPU; __SSE2__ is defined without any -m flag)
//   * AArch64 -> NEON (guaranteed on every ARMv8 / Apple Silicon CPU)
//   * anything else -> a plain scalar loop
// Each kernel produces bit-for-bit the same result as its scalar version, so callers stay byte-identical
// to VarDict regardless of which backend is compiled in.
#include <cstddef>
#include <cctype>

#if defined(__SSE2__)
  #include <emmintrin.h>
  #define VDCPP_SIMD_SSE2 1
#elif defined(__ARM_NEON) || defined(__aarch64__)
  #include <arm_neon.h>
  #define VDCPP_SIMD_NEON 1
#endif

namespace vardict {
namespace simd {

// Name of the compiled-in backend, for --version / diagnostics.
inline const char* backend() {
#if defined(VDCPP_SIMD_SSE2)
    return "SSE2";
#elif defined(VDCPP_SIMD_NEON)
    return "NEON";
#else
    return "scalar";
#endif
}

// ASCII upper-case in place: every byte in ['a','z'] has 0x20 cleared, all others untouched. This is
// exactly what `for (c : s) c = toupper((unsigned char)c)` does for ASCII FASTA bases (the only bytes
// faidx returns), 16 bytes at a time.
inline void toupper_ascii(char* p, size_t n) {
    size_t i = 0;
#if defined(VDCPP_SIMD_SSE2)
    const __m128i off = _mm_set1_epi8((char)0x80);          // map to signed space for unsigned compare
    const __m128i lo  = _mm_set1_epi8((char)('a' - 1 - 0x80));
    const __m128i hi  = _mm_set1_epi8((char)('z' + 1 - 0x80));
    const __m128i bit = _mm_set1_epi8((char)0x20);
    for (; i + 16 <= n; i += 16) {
        __m128i v  = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p + i));
        __m128i vo = _mm_add_epi8(v, off);
        __m128i ge = _mm_cmpgt_epi8(vo, lo);                // v >= 'a'
        __m128i le = _mm_cmpgt_epi8(hi, vo);                // v <= 'z'
        __m128i m  = _mm_and_si128(ge, le);
        v = _mm_sub_epi8(v, _mm_and_si128(m, bit));         // clear 0x20 where in range
        _mm_storeu_si128(reinterpret_cast<__m128i*>(p + i), v);
    }
#elif defined(VDCPP_SIMD_NEON)
    const uint8x16_t lo  = vdupq_n_u8('a');
    const uint8x16_t hi  = vdupq_n_u8('z');
    const uint8x16_t bit = vdupq_n_u8(0x20);
    for (; i + 16 <= n; i += 16) {
        uint8x16_t v = vld1q_u8(reinterpret_cast<const uint8_t*>(p + i));
        uint8x16_t m = vandq_u8(vcgeq_u8(v, lo), vcleq_u8(v, hi));
        v = vsubq_u8(v, vandq_u8(m, bit));
        vst1q_u8(reinterpret_cast<uint8_t*>(p + i), v);
    }
#endif
    for (; i < n; ++i) {
        unsigned char c = (unsigned char)p[i];
        if (c >= 'a' && c <= 'z') p[i] = (char)(c - 0x20);
    }
}

// Classify each base to its VarDict seed code: A->0, C->1, G->2, T->3, anything else (incl. N)->4.
// out must have room for n bytes. Bases are compared upper-case (the reference is upper-cased on load).
// Identical result to a per-byte lookup; used to feed the rolling k-mer encoder in the seed index.
inline void encode_bases(const char* s, unsigned char* out, size_t n) {
    size_t i = 0;
#if defined(VDCPP_SIMD_SSE2)
    const __m128i cA = _mm_set1_epi8('A'), cC = _mm_set1_epi8('C');
    const __m128i cG = _mm_set1_epi8('G'), cT = _mm_set1_epi8('T');
    const __m128i v0 = _mm_setzero_si128(), v1 = _mm_set1_epi8(1);
    const __m128i v2 = _mm_set1_epi8(2), v3 = _mm_set1_epi8(3), v4 = _mm_set1_epi8(4);
    for (; i + 16 <= n; i += 16) {
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(s + i));
        __m128i eA = _mm_cmpeq_epi8(v, cA), eC = _mm_cmpeq_epi8(v, cC);
        __m128i eG = _mm_cmpeq_epi8(v, cG), eT = _mm_cmpeq_epi8(v, cT);
        __m128i r = v4;                                             // default 4
        r = _mm_or_si128(_mm_and_si128(eA, v0), _mm_andnot_si128(eA, r));
        r = _mm_or_si128(_mm_and_si128(eC, v1), _mm_andnot_si128(eC, r));
        r = _mm_or_si128(_mm_and_si128(eG, v2), _mm_andnot_si128(eG, r));
        r = _mm_or_si128(_mm_and_si128(eT, v3), _mm_andnot_si128(eT, r));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(out + i), r);
    }
#elif defined(VDCPP_SIMD_NEON)
    const uint8x16_t cA = vdupq_n_u8('A'), cC = vdupq_n_u8('C');
    const uint8x16_t cG = vdupq_n_u8('G'), cT = vdupq_n_u8('T');
    for (; i + 16 <= n; i += 16) {
        uint8x16_t v = vld1q_u8(reinterpret_cast<const uint8_t*>(s + i));
        uint8x16_t r = vdupq_n_u8(4);
        r = vbslq_u8(vceqq_u8(v, cA), vdupq_n_u8(0), r);
        r = vbslq_u8(vceqq_u8(v, cC), vdupq_n_u8(1), r);
        r = vbslq_u8(vceqq_u8(v, cG), vdupq_n_u8(2), r);
        r = vbslq_u8(vceqq_u8(v, cT), vdupq_n_u8(3), r);
        vst1q_u8(out + i, r);
    }
#endif
    for (; i < n; ++i) {
        char ch = s[i];
        out[i] = ch == 'A' ? 0 : ch == 'C' ? 1 : ch == 'G' ? 2 : ch == 'T' ? 3 : 4;
    }
}

} // namespace simd
} // namespace vardict
