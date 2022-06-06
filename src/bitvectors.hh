/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef BITVECTOR_H
#define BITVECTOR_H

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

namespace bit
{
    inline __int512_t add_512bit(__int512_t a, __int512_t b)
    {
        __int512_t sum;
        char carry = 0;
        for (int i = 0; i < 8; i++)
            carry =
                _addcarry_u64(carry, a.words[i], b.words[i], &(sum.words[i]));
        return sum;
    }

    inline __int576_t add_576bit(__int576_t a, __int576_t b)
    {
        __int576_t sum;
        char carry = 0;
        for (int i = 0; i < 9; i++)
            carry =
                _addcarry_u64(carry, a.words[i], b.words[i], &(sum.words[i]));
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

    inline __int256_t mul_128bit(__int128_t a, __int128_t b)
    {
        __int128_t mask64b = 0xFFFFFFFFFFFFFFFFull;

        __int128_t ah = a >> 64;
        __int128_t al = a & mask64b;
        __int128_t bh = b >> 64;
        __int128_t bl = b & mask64b;

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
}
#endif
