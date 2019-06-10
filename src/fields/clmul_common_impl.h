/**********************************************************************
 * Copyright (c) 2018 Pieter Wuille, Greg Maxwell, Gleb Naumenko      *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _MINISKETCH_FIELDS_CLMUL_COMMON_IMPL_H_
#define _MINISKETCH_FIELDS_CLMUL_COMMON_IMPL_H_ 1

#include <stdint.h>
#include <x86intrin.h>

#include "../int_utils.h"
#include "../lintrans.h"

namespace {

template<typename I, int BITS, I MOD> I MulWithClMulReduce(I a, I b)
{
    static constexpr I MASK = Mask<BITS, I>();

    const __m128i MOD128 = _mm_cvtsi64_si128(MOD);
    __m128i product = _mm_clmulepi64_si128(_mm_cvtsi64_si128((uint64_t)a), _mm_cvtsi64_si128((uint64_t)b), 0x00);
    if (BITS <= 32) {
        __m128i high1 = _mm_srli_epi64(product, BITS);
        __m128i red1 = _mm_clmulepi64_si128(high1, MOD128, 0x00); 
        __m128i high2 = _mm_srli_epi64(red1, BITS);
        __m128i red2 = _mm_clmulepi64_si128(high2, MOD128, 0x00);
        return _mm_cvtsi128_si64(_mm_xor_si128(_mm_xor_si128(product, red1), red2)) & MASK;
    } else if (BITS == 64) {
        __m128i red1 = _mm_clmulepi64_si128(product, MOD128, 0x01);
        __m128i red2 = _mm_clmulepi64_si128(red1, MOD128, 0x01);
        return _mm_cvtsi128_si64(_mm_xor_si128(_mm_xor_si128(product, red1), red2));
    } else if ((BITS % 8) == 0) {
        __m128i high1 = _mm_srli_si128(product, BITS / 8);
        __m128i red1 = _mm_clmulepi64_si128(high1, MOD128, 0x00);
        __m128i high2 = _mm_srli_si128(red1, BITS / 8);
        __m128i red2 = _mm_clmulepi64_si128(high2, MOD128, 0x00);
        return _mm_cvtsi128_si64(_mm_xor_si128(_mm_xor_si128(product, red1), red2)) & MASK;
    } else {
        __m128i high1 = _mm_or_si128(_mm_srli_epi64(product, BITS), _mm_srli_si128(_mm_slli_epi64(product, 64 - BITS), 8));
        __m128i red1 = _mm_clmulepi64_si128(high1, MOD128, 0x00);
        if ((uint64_t(MOD) >> (66 - BITS)) == 0) {
            __m128i high2 = _mm_srli_epi64(red1, BITS);
            __m128i red2 = _mm_clmulepi64_si128(high2, MOD128, 0x00);
            return _mm_cvtsi128_si64(_mm_xor_si128(_mm_xor_si128(product, red1), red2)) & MASK;
        } else {
            __m128i high2 = _mm_or_si128(_mm_srli_epi64(red1, BITS), _mm_srli_si128(_mm_slli_epi64(red1, 64 - BITS), 8));
            __m128i red2 = _mm_clmulepi64_si128(high2, MOD128, 0x00);
            return _mm_cvtsi128_si64(_mm_xor_si128(_mm_xor_si128(product, red1), red2)) & MASK;
        }
    }
}

template<typename I, int BITS, int POS> I MulTrinomial(I a, I b)
{
    static constexpr I MASK = Mask<BITS, I>();

    __m128i product = _mm_clmulepi64_si128(_mm_cvtsi64_si128((uint64_t)a), _mm_cvtsi64_si128((uint64_t)b), 0x00);
    if (BITS <= 32) {
        __m128i high1 = _mm_srli_epi64(product, BITS);
        __m128i red1 = _mm_xor_si128(high1, _mm_slli_epi64(high1, POS));
        if (POS == 1) {
            return _mm_cvtsi128_si64(_mm_xor_si128(product, red1)) & MASK;
        } else {
            __m128i high2 = _mm_srli_epi64(red1, BITS);
            __m128i red2 = _mm_xor_si128(high2, _mm_slli_epi64(high2, POS));
            return _mm_cvtsi128_si64(_mm_xor_si128(_mm_xor_si128(product, red1), red2)) & MASK;
        }
    } else {
        __m128i high1 = _mm_or_si128(_mm_srli_epi64(product, BITS), _mm_srli_si128(_mm_slli_epi64(product, 64 - BITS), 8));
        if (BITS + POS <= 66) {
            __m128i red1 = _mm_xor_si128(high1, _mm_slli_epi64(high1, POS));
            if (POS == 1) {
                return _mm_cvtsi128_si64(_mm_xor_si128(product, red1)) & MASK;
            } else if (BITS + POS <= 66) {
                __m128i high2 = _mm_srli_epi64(red1, BITS);
                __m128i red2 = _mm_xor_si128(high2, _mm_slli_epi64(high2, POS));
                return _mm_cvtsi128_si64(_mm_xor_si128(_mm_xor_si128(product, red1), red2)) & MASK;
            }
        } else {
            const __m128i MOD128 = _mm_cvtsi64_si128(1 + (((uint64_t)1) << POS));
            __m128i red1 = _mm_clmulepi64_si128(high1, MOD128, 0x00);
            __m128i high2 = _mm_or_si128(_mm_srli_epi64(red1, BITS), _mm_srli_si128(_mm_slli_epi64(red1, 64 - BITS), 8));
            __m128i red2 = _mm_xor_si128(high2, _mm_slli_epi64(high2, POS));
            return _mm_cvtsi128_si64(_mm_xor_si128(_mm_xor_si128(product, red1), red2)) & MASK;
        }
    }
}

/** Implementation of fields that use the SSE clmul intrinsic for multiplication. */
template<typename I, int B, I MOD, I (*MUL)(I, I), typename F, const F* SQR, const F* SQR2, const F* SQR4, const F* SQR8, const F* SQR16, const F* QRT, typename T, const T* LOAD, const T* SAVE> struct GenField
{
    typedef BitsInt<I, B> O;
    typedef LFSR<O, MOD> L;

    I m_val;
    explicit constexpr GenField(I val) : m_val(val) {}

    static inline I Mul(I a, I b) { return MUL(a, b); }

    static inline constexpr I Sqr1(I a) { return SQR->template Map<O>(a); }
    static inline constexpr I Sqr2(I a) { return SQR2->template Map<O>(a); }
    static inline constexpr I Sqr4(I a) { return SQR4->template Map<O>(a); }
    static inline constexpr I Sqr8(I a) { return SQR8->template Map<O>(a); }
    static inline constexpr I Sqr16(I a) { return SQR16->template Map<O>(a); }

public:
    static constexpr int BITS = B;

    static inline constexpr GenField Zero() { return GenField(0); }
    static inline constexpr GenField One() { return GenField(1); }

    inline constexpr GenField() : m_val(0) {}
    inline constexpr bool IsZero() const { return m_val == 0; }
    inline constexpr bool IsOne() const { return m_val == 1; }

    inline constexpr bool friend operator==(GenField a, GenField b) { return a.m_val == b.m_val; }
    inline constexpr bool friend operator!=(GenField a, GenField b) { return a.m_val != b.m_val; }
    inline constexpr bool friend operator<(GenField a, GenField b) { return a.m_val < b.m_val; }

    inline friend constexpr GenField operator+(GenField a, GenField b) { return GenField(a.m_val ^ b.m_val); }

    inline GenField& operator+=(GenField a) { m_val ^= a.m_val; return *this; }

    inline constexpr GenField Mul2() const { return GenField(L::Call(m_val)); }

    inline friend GenField operator*(GenField a, GenField b) { return GenField(Mul(a.m_val, b.m_val)); }

    class Multiplier
    {
        I m_val;
    public:
        inline constexpr explicit Multiplier(GenField a) : m_val(a.m_val) {}
        constexpr GenField operator()(GenField a) const { return GenField(Mul(m_val, a.m_val)); }
    };

    /** Compute the square of a. */
    inline constexpr GenField Sqr() const { return GenField(Sqr1(m_val)); }

    /** Compute x such that x^2 + x = a (undefined result if no solution exists). */
    inline constexpr GenField Qrt() const { return GenField(QRT->template Map<O>(m_val)); }

    /** Compute the inverse of x1. */
    inline GenField Inv() const { return GenField(InvLadder<I, O, B, Mul, Sqr1, Sqr2, Sqr4, Sqr8, Sqr16>(m_val)); }

    /** Generate a random field element. */
    static GenField FromSeed(uint64_t seed) {
        uint64_t k0 = 0x434c4d554c466c64ull; // "CLMULFld"
        uint64_t k1 = seed;
        uint64_t count = ((uint64_t)BITS) << 32;
        I ret;
        do {
            ret = O::Mask(I(SipHash(k0, k1, count++)));
        } while(ret == 0);
        return GenField(LOAD->template Map<O>(ret));
    }

    static GenField Deserialize(BitReader& in) { return GenField(LOAD->template Map<O>(in.Read<B, I>())); }

    void Serialize(BitWriter& out) const { out.Write<B, I>(SAVE->template Map<O>(m_val)); }

    static constexpr GenField FromUint64(uint64_t x) { return GenField(LOAD->template Map<O>(O::Mask(I(x)))); }
    constexpr uint64_t ToUint64() const { return uint64_t(SAVE->template Map<O>(m_val)); }
};

template<typename I, int B, I MOD, typename F, const F* SQR, const F* SQR2, const F* SQR4, const F* SQR8, const F* SQR16, const F* QRT, typename T, const T* LOAD, const T* SAVE>
using Field = GenField<I, B, MOD, MulWithClMulReduce<I, B, MOD>, F, SQR, SQR2, SQR4, SQR8, SQR16, QRT, T, LOAD, SAVE>;

template<typename I, int B, int POS, typename F, const F* SQR, const F* SQR2, const F* SQR4, const F* SQR8, const F* SQR16, const F* QRT, typename T, const T* LOAD, const T* SAVE>
using FieldTri = GenField<I, B, I(1) + (I(1) << POS), MulTrinomial<I, B, POS>, F, SQR, SQR2, SQR4, SQR8, SQR16, QRT, T, LOAD, SAVE>;

}

#endif
