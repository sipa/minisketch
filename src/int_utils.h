/**********************************************************************
 * Copyright (c) 2018 Pieter Wuille, Greg Maxwell, Gleb Naumenko      *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _MINISKETCH_INT_UTILS_H_
#define _MINISKETCH_INT_UTILS_H_

#include <stdlib.h>

#include <limits>
#include <algorithm>
#include <type_traits>

template<int bits>
static constexpr inline uint64_t Rot(uint64_t x) { return (x << bits) | (x >> (64 - bits)); }

static inline void SipHashRound(uint64_t& v0, uint64_t& v1, uint64_t& v2, uint64_t& v3) {
    v0 += v1; v1 = Rot<13>(v1); v1 ^= v0;
    v0 = Rot<32>(v0);
    v2 += v3; v3 = Rot<16>(v3); v3 ^= v2;
    v0 += v3; v3 = Rot<21>(v3); v3 ^= v0;
    v2 += v1; v1 = Rot<17>(v1); v1 ^= v2;
    v2 = Rot<32>(v2);
}

inline uint64_t SipHash(uint64_t k0, uint64_t k1, uint64_t data) {
    uint64_t v0 = 0x736f6d6570736575ULL ^ k0;
    uint64_t v1 = 0x646f72616e646f6dULL ^ k1;
    uint64_t v2 = 0x6c7967656e657261ULL ^ k0;
    uint64_t v3 = 0x7465646279746573ULL ^ k1 ^ data;
    SipHashRound(v0, v1, v2, v3);
    SipHashRound(v0, v1, v2, v3);
    v0 ^= data;
    v3 ^= 0x800000000000000ULL;
    SipHashRound(v0, v1, v2, v3);
    SipHashRound(v0, v1, v2, v3);
    v0 ^= 0x800000000000000ULL;
    v2 ^= 0xFF;
    SipHashRound(v0, v1, v2, v3);
    SipHashRound(v0, v1, v2, v3);
    SipHashRound(v0, v1, v2, v3);
    SipHashRound(v0, v1, v2, v3);
    return v0 ^ v1 ^ v2 ^ v3;
}

class BitWriter {
    unsigned char state = 0;
    int offset = 0;
    unsigned char* out;

public:
    BitWriter(unsigned char* output) : out(output) {}

    template<int BITS, typename I>
    inline void Write(I val) {
        int bits = BITS;
        if (bits + offset >= 8) {
            state |= ((val & ((I(1) << (8 - offset)) - 1)) << offset);
            *(out++) = state;
            val >>= (8 - offset);
            bits -= 8 - offset;
            offset = 0;
            state = 0;
        }
        while (bits >= 8) {
            *(out++) = val & 255;
            val >>= 8;
            bits -= 8;
        }
        state |= ((val & ((I(1) << bits) - 1)) << offset);
        offset += bits;
    }

    inline void Flush() {
        if (offset) {
            *(out++) = state;
            state = 0;
            offset = 0;
        }
    }
};

class BitReader {
    unsigned char state = 0;
    int offset = 0;
    const unsigned char* in;

public:
    BitReader(const unsigned char* input) : in(input) {}

    template<int BITS, typename I>
    inline I Read() {
        int bits = BITS;
        if (offset >= bits) {
            I ret = state & ((1 << bits) - 1);
            state >>= bits;
            offset -= bits;
            return ret;
        }
        I val = state;
        int out = offset;
        while (out + 8 <= bits) {
            val |= ((I(*(in++))) << out);
            out += 8;
        }
        if (out < bits) {
            unsigned char c = *(in++);
            val |= (c & ((I(1) << (bits - out)) - 1)) << out;
            state = c >> (bits - out);
            offset = 8 - (bits - out);
        } else {
            state = 0;
            offset = 0;
        }
        return val;
    }
};

/** Return a value of type I with its `bits` lowest bits set (bits must be > 0). */
template<int BITS, typename I>
constexpr inline I Mask() { return ((I((I(-1)) << (std::numeric_limits<I>::digits - BITS))) >> (std::numeric_limits<I>::digits - BITS)); }

template<typename I>
struct DoubleInt
{
    I low, high;

    explicit DoubleInt(I val) : low(val), high(0) {}
};

template<typename I, int BITS>
class BitsInt {
private:
    static_assert(std::is_integral<I>::value && std::is_unsigned<I>::value, "BitsInt requires an unsigned integer type");
    static_assert(BITS > 0 && BITS <= std::numeric_limits<I>::digits, "BitsInt requires 1 <= Bits <= representation type size");

    static constexpr I MASK = Mask<BITS, I>();

public:

    typedef I Repr;

    static constexpr int SIZE = BITS;

    static void inline Swap(I& a, I& b) {
        std::swap(a, b);
    }

    static constexpr inline bool IsZero(I a) { return a == 0; }
    static constexpr inline I Mask(I val) { return val & MASK; }
    static constexpr inline I Shift(I val, int bits) { return ((val << bits) & MASK); }
    static constexpr inline I UnsafeShift(I val, int bits) { return (val << bits); }

    template<int Offset, int Count>
    static constexpr inline int MidBits(I val) {
        static_assert(Count > 0, "BITSInt::MidBits needs Count > 0");
        static_assert(Count + Offset <= BITS, "BitsInt::MidBits overflow of Count+Offset");
        return (val >> Offset) & ((I(1) << Count) - 1);
    }

    template<int Count>
    static constexpr inline int TopBits(I val) {
        static_assert(Count > 0, "BitsInt::TopBits needs Count > 0");
        static_assert(Count <= BITS, "BitsInt::TopBits needs Offset <= BITS");
        return val >> (BITS - Count);
    }

    static inline constexpr I CondXorWith(I val, bool cond, I v) {
        return val ^ (-I(cond) & v);
    }

    template<I MOD>
    static inline constexpr I CondXorWith(I val, bool cond) {
        return val ^ (-I(cond) & MOD);
    }

    static inline int Bits(I val, int max) {
#ifdef HAVE_CLZ
        (void)max;
        if (val == 0) return 0;
        if (std::numeric_limits<unsigned>::digits >= std::numeric_limits<I>::digits) {
            return std::numeric_limits<unsigned>::digits - __builtin_clz(val);
        } else if (std::numeric_limits<unsigned long>::digits >= std::numeric_limits<I>::digits) {
            return std::numeric_limits<unsigned long>::digits - __builtin_clzl(val);
        } else {
            return std::numeric_limits<unsigned long long>::digits - __builtin_clzll(val);
        }
#else
        while (max && (val >> (max - 1) == 0)) --max;
        return max;
#endif
    }
};

template<typename I, int LBITS, int HBITS>
class BitsDoubleInt {
private:
    static_assert(std::is_integral<I>::value && std::is_unsigned<I>::value, "DoubleBitsInt requires an unsigned integer type");
    static_assert(LBITS > 0 && HBITS >0 && LBITS <= std::numeric_limits<I>::digits && HBITS <= std::numeric_limits<I>::digits, "DoubleBits requires 1 <= {LBITS, HBITS} <= representation type size");

    static constexpr I LMASK = Mask<LBITS, I>();
    static constexpr I HMASK = Mask<HBITS, I>();

public:

    typedef DoubleInt<I> Repr;

    static constexpr int SIZE = LBITS + HBITS;

    static void inline Swap(Repr& a, Repr& b) {
        std::swap(a.low, b.low);
        std::swap(a.high, b.high);
    }

    static constexpr inline bool IsZero(const Repr& a) { return a.low == 0 && a.high == 0; }
    static constexpr inline Repr Mask(const Repr& a) { return Repr{a.low & LMASK, a.high & HMASK}; }
    static constexpr inline Repr Shift(const Repr& a, int bits) { return Repr{(a.low << bits) & LMASK, ((a.high << bits) | (a.low >> (LBITS - bits))) & HMASK}; }
    static constexpr inline Repr UnsafeShift(const Repr& a, int bits) { return Repr{(a.low << bits) & LMASK, (a.high << bits) | (a.low >> (LBITS - bits))}; }

    template<int Offset, int Count>
    static constexpr inline int MidBits(const Repr& val) {
        static_assert(Count > 0, "Bla");
        static_assert(Count + Offset <= LBITS + HBITS, "Blu");
        return ((Count + Offset <= LBITS) ?
            (val.low >> Offset) :
            ((Offset >= LBITS) ?
                (val.high >> (Offset - LBITS)) :
                ((val.low >> Offset) | (val.high << (LBITS - Offset))))) & ((I(1) << Count) - 1);
    }

    template<int Count>
    static constexpr inline int TopBits(const Repr& val) {
        static_assert(Count > 0, "Bli");
        static_assert(Count <= HBITS, "Ble");
        return (val.high >> (HBITS - Count));
    }

    static inline constexpr Repr CondXorWith(const Repr& val, bool cond, const Repr& v) {
        return Repr{val.low ^ (-I(cond) & v.low), val.high ^ (-I(cond) & v.high)};
    }

    template<I MOD>
    static inline constexpr I CondXorWith(const Repr& val, bool cond) {
        return Repr{val.low ^ (-I(cond) & MOD), val.high};
    }

    static inline int Bits(const Repr& val, int max) {
        if (max > LBITS && val.high) {
            return BitsInt<I, HBITS>::Bits(val.high, max - LBITS) + LBITS;
        }
        return BitsInt<I, LBITS>::Bits(val.low, max);
    }
};

/** Class which implements a stateless LFSR for generic moduli. */
template<typename F, uint32_t MOD>
struct LFSR {
    typedef typename F::Repr I;
    /** Shift a value `a` up once, treating it as an `N`-bit LFSR, with pattern `MOD`. */
    static inline constexpr I Call(const I& a) {
        return F::template CondXorWith<MOD>(F::Shift(a, 1), F::template TopBits<1>(a));
    }
};

/** Helper class for carryless multiplications. */
template<typename I, int N, typename L, typename F, int K> struct GFMulHelper;
template<typename I, int N, typename L, typename F> struct GFMulHelper<I, N, L, F, 0>
{
    static inline constexpr I Run(const I& a, const I& b) { return I(); }
};
template<typename I, int N, typename L, typename F, int K> struct GFMulHelper
{
    static inline constexpr I Run(const I& a, const I& b) { return F::CondXorWith(GFMulHelper<I, N, L, F, K - 1>::Run(L::Call(a), b), F::template MidBits<N - K, 1>(b), a); }
};

/** Compute the carry-less multiplication of a and b, with N bits, using L as LFSR type. */
template<typename I, int N, typename L, typename F> inline constexpr I GFMul(const I& a, const I& b) { return GFMulHelper<I, N, L, F, N>::Run(a, b); }

/** Compute the inverse of x using an extgcd algorithm. */
template<typename I, typename F, int BITS, uint32_t MOD>
inline I InvExtGCD(I x)
{
    if (F::IsZero(x)) return x;
    I t(0), newt(1);
    I r(MOD), newr = x;
    int rlen = BITS + 1, newrlen = F::Bits(newr, BITS);
    while (!F::IsZero(newr)) {
        int q = rlen - newrlen;
        r ^= F::Shift(newr, q);
        t ^= F::UnsafeShift(newt, q);
        rlen = F::Bits(r, rlen - 1);
        if (r < newr) {
            F::Swap(t, newt);
            F::Swap(r, newr);
            std::swap(rlen, newrlen);
        }
    }
    return t;
}

/** Compute the inverse of x1 using an exponentiation ladder.
 *
 * The `MUL` argument is a multiplication function, `SQR` is a squaring function, and the `SQRi` arguments
 * compute x**(2**i).
 */
template<typename I, typename F, int BITS, I (*MUL)(I, I), I (*SQR)(I), I (*SQR2)(I), I(*SQR4)(I), I(*SQR8)(I), I(*SQR16)(I)>
inline I InvLadder(I x1)
{
    static constexpr int INV_EXP = BITS - 1;
    I x2 = (INV_EXP >= 2) ? MUL(SQR(x1), x1) : I(0);
    I x4 = (INV_EXP >= 4) ? MUL(SQR2(x2), x2) : I(0);
    I x8 = (INV_EXP >= 8) ? MUL(SQR4(x4), x4) : I(0);
    I x16 = (INV_EXP >= 16) ? MUL(SQR8(x8), x8) : I(0);
    I x32 = (INV_EXP >= 32) ? MUL(SQR16(x16), x16) : I(0);
    I r;
    if (INV_EXP >= 32) {
        r = x32;
    } else if (INV_EXP >= 16) {
        r = x16;
    } else if (INV_EXP >= 8) {
        r = x8;
    } else if (INV_EXP >= 4) {
        r = x4;
    } else if (INV_EXP >= 2) {
        r = x2;
    } else {
        r = x1;
    }
    if (INV_EXP >= 32 && (INV_EXP & 16)) r = MUL(SQR16(r), x16);
    if (INV_EXP >= 16 && (INV_EXP & 8)) r = MUL(SQR8(r), x8);
    if (INV_EXP >= 8 && (INV_EXP & 4)) r = MUL(SQR4(r), x4);
    if (INV_EXP >= 4 && (INV_EXP & 2)) r = MUL(SQR2(r), x2);
    if (INV_EXP >= 2 && (INV_EXP & 1)) r = MUL(SQR(r), x1);
    return SQR(r);
}

#endif
