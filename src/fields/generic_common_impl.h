/**********************************************************************
 * Copyright (c) 2018 Pieter Wuille, Greg Maxwell, Gleb Naumenko      *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _MINISKETCH_FIELDS_GENERIC_COMMON_IMPL_H_
#define _MINISKETCH_FIELDS_GENERIC_COMMON_IMPL_H_ 1

#include <stdint.h>
#include <random>

#include "../int_utils.h"
#include "../lintrans.h"

namespace {

/** Generic implementation for fields whose elements can be represented by an integer type. */
template<typename I, int B, uint32_t MOD, typename F, typename T, const F* SQR, const F* QRT> class Field
{
    typedef BitsInt<I, B> O;
    typedef LFSR<O, MOD> L;

    I m_val;
    explicit constexpr Field(I val) : m_val(val) {}

public:
    static constexpr int BITS = B;

    static constexpr Field Zero() { return Field(I(0)); }
    static constexpr Field One() { return Field(I(1)); }

    constexpr inline Field() : m_val(0) {}
    constexpr inline bool IsZero() const { return m_val == 0; }
    constexpr inline bool IsOne() const { return m_val == 1; }

    constexpr inline bool friend operator==(Field a, Field b) { return a.m_val == b.m_val; }
    constexpr inline bool friend operator!=(Field a, Field b) { return a.m_val != b.m_val; }
    constexpr inline bool friend operator<(Field a, Field b) { return a.m_val < b.m_val; }

    constexpr inline Field Mul2() const { return Field(L::Call(m_val)); }

    class Multiplier
    {
        T table;
    public:
        explicit Multiplier(Field a) { table.template Build<L::Call>(a.m_val); }
        constexpr inline Field operator()(Field a) const { return Field(table.template Map<O>(a.m_val)); }
    };

    inline friend constexpr Field operator+(Field a, Field b) { return Field(a.m_val ^ b.m_val); }
    inline Field& operator+=(Field a) { m_val ^= a.m_val; return *this; }
    friend Field operator*(Field a, Field b) { return Field(GFMul<I, B, L, O>(a.m_val, b.m_val)); }

    /** Compute the square of a. */
    inline constexpr Field Sqr() const { return Field(SQR->template Map<O>(m_val)); }

    /** Compute x such that x^2 + x = a (undefined result if no solution exists). */
    inline constexpr Field Qrt() const { return Field(QRT->template Map<O>(m_val)); }

    /** Compute the inverse of x1. */
    Field Inv() const { return Field(InvExtGCD<I, O, B, MOD>(m_val)); }

    /** Generate a random field element. */
    static Field FromSeed(uint64_t seed) {
        uint64_t k0 = 0x496e744669656c64ull; // "IntField"
        uint64_t k1 = seed;
        uint64_t count = ((uint64_t)BITS) << 32;
        I ret;
        do {
            ret = O::Mask(I(SipHash(k0, k1, count++)));
        } while(ret == 0);
        return Field(ret);
    }

    static Field Deserialize(BitReader& in) { return Field(in.template Read<BITS, I>()); }

    void Serialize(BitWriter& out) const { out.template Write<BITS, I>(m_val); }

    static constexpr Field FromUint64(uint64_t x) { return Field(O::Mask(I(x))); }
    constexpr uint64_t ToUint64() const { return uint64_t(m_val); }
};

/** Generic implementation for fields whose elements can be represented by an integer type. */
template<typename II, int LB, int HB, uint32_t MOD, typename F, typename T, const F* SQR, const F* QRT> class DoubleField
{
    typedef DoubleInt<II> I;
    typedef BitsDoubleInt<II, LB, HB> O;
    typedef LFSR<O, MOD> L;

    I m_val;
    explicit constexpr Field(I val) : m_val(val) {}

public:
    static constexpr int BITS = LB + HB;

    static constexpr Field Zero() { return Field(I(0)); }
    static constexpr Field One() { return Field(I(1)); }

    constexpr inline Field() : m_val(0) {}
    constexpr inline bool IsZero() const { return m_val.low == 0 && m_val.high == 0; }
    constexpr inline bool IsOne() const { return m_val.low == 1 && m_val.high == 0; }

    constexpr inline bool friend operator==(Field a, Field b) { return a.m_val.low == b.m_val.low && a.m_val.high == b.m_val.high; }
    constexpr inline bool friend operator!=(Field a, Field b) { return a.m_val.low != b.m_val.low || a.m_val.high != b.m_val.high; }
    constexpr inline bool friend operator<(Field a, Field b) { return a.m_val.high < b.m_val.high || (a.m_val.high == b.m_val.high && a.m_val.low < b.m_val.low); }

    constexpr inline Field Mul2() const { return Field(L::Call(m_val)); }

    class Multiplier
    {
        T table;
    public:
        explicit Multiplier(Field a) { table.template Build<L::Call>(a.m_val); }
        constexpr inline Field operator()(Field a) const { return Field(table.template Map<O>(a.m_val)); }
    };

    inline friend constexpr Field operator+(Field a, Field b) { return Field(a.m_val ^ b.m_val); }
    inline Field& operator+=(Field a) { m_val ^= a.m_val; return *this; }
    friend Field operator*(Field a, Field b) { return Field(GFMul<I, B, L, O>(a.m_val, b.m_val)); }

    /** Compute the square of a. */
    inline constexpr Field Sqr() const { return Field(SQR->template Map<O>(m_val)); }

    /** Compute x such that x^2 + x = a (undefined result if no solution exists). */
    inline constexpr Field Qrt() const { return Field(QRT->template Map<O>(m_val)); }

    /** Compute the inverse of x1. */
    Field Inv() const { return Field(InvExtGCD<I, O, B, MOD>(m_val)); }

    /** Generate a random field element. */
    static Field FromSeed(uint64_t seed) {
        uint64_t k0 = 0x496e744669656c64ull; // "IntField"
        uint64_t k1 = seed;
        uint64_t count = ((uint64_t)BITS) << 32;
        I ret;
        do {
            ret = O::Mask(I(SipHash(k0, k1, count++)));
        } while(ret == 0);
        return Field(ret);
    }

    static Field Deserialize(BitReader& in) { return Field(in.template Read<BITS, I>()); }

    void Serialize(BitWriter& out) const { out.template Write<BITS, I>(m_val); }

    static constexpr Field FromUint64(uint64_t x) { return Field(O::Mask(I(x))); }
    constexpr uint64_t ToUint64() const { return uint64_t(m_val); }
};

}

#endif
