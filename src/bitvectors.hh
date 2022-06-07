/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef BITVECTOR_H
#define BITVECTOR_H

struct uint128_t
{
    uint64_t words[2];
};

struct __int256_t
{
    long long unsigned int words[4];
};

struct __int512_t
{
    long long unsigned int words[8];
};

struct __int576_t
{
    long long unsigned int words[9];
};

#define ADD_WORDS(n)                                                           \
{                                                                              \
    char carry = 0;                                                            \
    for (int i = 0; i < n; i++)                                                \
        carry = _addcarry_u64(carry, a.words[i], b.words[i], &(sum.words[i])); \
}

namespace bit
{
    inline __int512_t add_512bit(__int512_t a, __int512_t b)
    {
        __int512_t sum;
        ADD_WORDS(8)
        return sum;
    }

    inline __int576_t pad_words(__int512_t a, int n)
    {
        __int576_t padded;
        for (int i = 0; i < n; i++)
            padded.words[i] = 0;

        for (int i = n; i < 9; i++)
            padded.words[i] = a.words[i - n];

        return padded;
    }

    inline __int576_t widen_512bits(__int512_t a)
    {
        __int576_t wide;
        for (int i = 0; i < 8; i++)
            wide.words[i] = a.words[i];

        wide.words[8] = 0;

        return wide;
    }

    inline __int576_t add_576bit(__int576_t a, __int576_t b)
    {
        __int576_t sum;
        ADD_WORDS(9)
        return sum;
    }

    /* pos <= 64 */
    inline __int576_t lshift_512bit(__int512_t a, int pos)
    {
        __int576_t shift;
        shift.words[0] = a.words[0] << pos;
        for (int i = 1; i < 8; i++)
        {
            shift.words[i] = a.words[i] << pos;
            shift.words[i] |= a.words[i-1] >> (64 - pos);
        }
        shift.words[8] = a.words[7] >> (64 - pos);
        return shift;
    }

    inline char add_128bit_carry(uint128_t a, uint128_t b, long long unsigned int *sum)
    {
        char carry = 0;
        carry = _addcarry_u64(carry, a.words[0], b.words[0], sum);
        carry = _addcarry_u64(carry, a.words[1], b.words[1], sum + 1);
        return carry;
    }

    inline __int256_t mul_128bit(uint128_t a, uint128_t b)
    {
        __int128_t mask64b = 0xFFFFFFFFFFFFFFFFull;

        __int128_t ah = a.words[1];
        __int128_t al = a.words[0];
        __int128_t bh = b.words[1];
        __int128_t bl = b.words[0];

        /* see https://stackoverflow.com/a/26855440 for logic. */

        __int128_t ahbh = ah*bh;
        __int128_t ahbl = ah*bl;
        __int128_t albh = al*bh;
        __int128_t albl = al*bl;

        __int128_t mid = (albl >> 64) + (ahbl & mask64b) + (albh & mask64b);

        __int256_t ret;
        ret.words[0] = albl & mask64b;
        ret.words[1] = mid & mask64b;
        __int128_t hi = ahbh + (ahbl >> 64) + (albh >> 64) + (mid >> 64);
        ret.words[2] = hi & mask64b;
        ret.words[3] = hi >> 64;

        return ret;
    }

    inline __int512_t mul_256bit_64bit(__int256_t a, uint64_t b)
    {
        __int256_t ahbl = bit::mul_128bit(
            { a.words[2], a.words[3] },
            { b, 0 }
        );

        __int256_t albl = bit::mul_128bit(
            { a.words[0], a.words[1] },
            { b, 0}
        );

        __int256_t mid;
        unsigned char carry = 0;
        long long unsigned sum[2];
        carry = add_128bit_carry(
            { albl.words[2], 0 },
            { ahbl.words[0], ahbl.words[1] },
            sum
        );

        mid.words[0] = sum[0];
        mid.words[1] = sum[1];
        mid.words[2] = carry;

        __int512_t ret;
        ret.words[0] = albl.words[0];
        ret.words[1] = albl.words[1];
        ret.words[2] = mid.words[0];
        ret.words[3] = mid.words[1];
        ret.words[4] = mid.words[2] + ahbl.words[2];
        ret.words[5] = 0;
        ret.words[6] = 0;
        ret.words[7] = 0;

        return ret;
    }

    inline __int512_t mul_256bit(__int256_t a, __int256_t b)
    {
        __int256_t ahbh = bit::mul_128bit(
            { a.words[2], a.words[3] },
            { b.words[2], b.words[3] }
        );
        __int256_t ahbl = bit::mul_128bit(
            { a.words[2], a.words[3] },
            { b.words[0], b.words[1] }
        );
        __int256_t albh = bit::mul_128bit(
            { a.words[0], a.words[1] },
            { b.words[2], b.words[3] }
        );
        __int256_t albl = bit::mul_128bit(
            { a.words[0], a.words[1] },
            { b.words[0], b.words[1] }
        );

        __int256_t mid;
        unsigned char carry = 0;
        long long unsigned int sum[2];
        carry = add_128bit_carry(
            { albl.words[2], albl.words[3] },
            { ahbl.words[0], ahbl.words[1] },
            sum
        );
        carry += add_128bit_carry(
            { sum[0], sum[1] },
            { albh.words[0], albh.words[1] },
            sum
        );
        mid.words[0] = sum[0];
        mid.words[1] = sum[1];
        mid.words[2] = carry;

        __int512_t ret;
        ret.words[0] = albl.words[0];
        ret.words[1] = albl.words[1];
        ret.words[2] = mid.words[0];
        ret.words[3] = mid.words[1];

        sum[0] = 0;
        sum[1] = 0;
        carry = add_128bit_carry(
            { ahbh.words[0], ahbh.words[1] },
            { ahbl.words[2], ahbl.words[3] },
            sum
        );
        carry += add_128bit_carry(
            { sum[0], sum[1] },
            { albh.words[2], albh.words[3] },
            sum
        );
        carry += add_128bit_carry(
            { sum[0], sum[1] },
            { mid.words[2], 0 },
            sum
        );

        ret.words[4] = sum[0];
        ret.words[5] = sum[1];

        sum[0] = 0;
        sum[1] = 0;
        carry = add_128bit_carry(
            { carry, 0 },
            { ahbh.words[2], ahbh.words[3] },
            sum
        );

        ret.words[6] = sum[0];
        ret.words[7] = sum[1];

        return ret;
    }
}
#endif
