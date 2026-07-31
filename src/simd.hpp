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

// Number of positions where a[i] != b[i] over the first n bytes (a Hamming-style mismatch count on
// equal-length sequences). Used by read/reference comparison hot paths; matches the scalar loop
// `for (i<n) cnt += a[i]!=b[i]` exactly.
inline int mismatch_count(const char* a, const char* b, size_t n) {
    size_t i = 0;
    int cnt = 0;
#if defined(VDCPP_SIMD_SSE2)
    for (; i + 16 <= n; i += 16) {
        __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i*>(a + i));
        __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i*>(b + i));
        // equal-byte mask -> bit per byte; popcount of the *cleared* bits = mismatches.
        unsigned eq = (unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(va, vb));
        cnt += 16 - __builtin_popcount(eq & 0xFFFFu);
    }
#elif defined(VDCPP_SIMD_NEON)
    for (; i + 16 <= n; i += 16) {
        uint8x16_t va = vld1q_u8(reinterpret_cast<const uint8_t*>(a + i));
        uint8x16_t vb = vld1q_u8(reinterpret_cast<const uint8_t*>(b + i));
        uint8x16_t ne = vmvnq_u8(vceqq_u8(va, vb));         // 0xFF where bytes differ
        // sum the 0xFF/0x00 lanes: (sum/255) = number of differing bytes.
        cnt += (int)(vaddvq_u8(ne) / 255u);
    }
#endif
    for (; i < n; ++i) cnt += (a[i] != b[i]);
    return cnt;
}

} // namespace simd
} // namespace vardict
