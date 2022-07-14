/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef P_FMATRIX_H
#define P_FMATRIX_H

#include <valarray>
#include <vector>
#include <immintrin.h>

#include "gf.hh"
#include "global.hh"
#include "fmatrix.hh"

#define VECTOR_N 8

#define DET_LOOP(index)                                         \
    {                                                           \
        int r0 = VECTOR_N*col + index;                          \
        __m256i mx = this->get(r0, col);                        \
        int mxi = r0;                                           \
        char cmpmsk = 1 << (2*(3-index));                       \
        for (int row = r0 + 1; row < this->rows; row++)         \
        {                                                       \
            char cmp = _mm256_cmp_epu32_mask(                   \
                mx,                                             \
                this->get(row, col),                            \
                0x1                                             \
            );                                                  \
            if (cmp & cmpmsk)                                   \
            {                                                   \
                mx = this->get(row,col);                        \
                mxi = row;                                      \
            }                                                   \
        }                                                       \
        uint64_t mx_ext =                                       \
            _mm256_extract_epi64(mx, 3 - index);                \
        if (mx_ext == 0)                                        \
            return global::F.zero();                            \
        if (mxi != r0)                                          \
            this->swap_rows(mxi, r0);                           \
        det = global::F.rem(                                    \
            global::F.clmul(det, mx_ext)                        \
        );                                                      \
        mx_ext = global::F.ext_euclid(mx_ext);                  \
        this->mul_row(r0, mx_ext);                              \
        for (int row = r0 + 1; row < this->rows; row++)         \
        {                                                       \
            /* use shuffle magix here? */                       \
            uint64_t val = _mm256_extract_epi64(                \
                this->get(row, col),                            \
                3 - index                                       \
            );                                                  \
            this->row_op(r0, row, val);                         \
        }                                                       \
    }

class Packed_FMatrix
{
private:
    int rows;
    int cols;
    // original matrix n moduloe VECTOR_N
    int nmod;
    std::vector<__m256i> m;

    __m256i get(int row, int col)
    {
        return this->m[row*this->cols + col];
    }

    void set(int row, int col, __m256i v)
    {
        this->m[row*this->cols + col] = v;
    }

public:
    Packed_FMatrix(
        const FMatrix &matrix
    )
    {
        this->nmod = matrix.get_n() % VECTOR_N;
        this->rows = matrix.get_n();
        if (this->rows % VECTOR_N)
            this->rows += VECTOR_N - (matrix.get_n() % VECTOR_N);
        this->cols = this->rows / VECTOR_N;

        this->m.resize(this->rows * this->cols);

        for (int r = 0; r < matrix.get_n(); r++)
        {
            for (int c = 0; c < matrix.get_n() / VECTOR_N; c++)
                this->set(r, c,
                          _mm256_set_epi64x(
                              matrix(r, VECTOR_N*c + 0).get_repr() << 32 | matrix(r, VECTOR_N*c + 1),
                              matrix(r, VECTOR_N*c + 2).get_repr() << 32 | matrix(r, VECTOR_N*c + 3),
                              matrix(r, VECTOR_N*c + 4).get_repr() << 32 | matrix(r, VECTOR_N*c + 5),
                              matrix(r, VECTOR_N*c + 6).get_repr() << 32 | matrix(r, VECTOR_N*c + 7)
                          )
                    );
            if (this->nmod)
            {
                int c = this->cols - 1;
                uint64_t elems[4];
                elems[0] = 0; elems[1] = 0; elems[2] = 0; elems[3] = 0;
                for (int i = 0; i < this->nmod; i++)
                    elems[i/2] |= matrix(r, VECTOR_N*c + i).get_repr() << (32*(1 - i%2));

                this->set(r, c, _mm256_set_epi64x(
                              elems[0],
                              elems[1],
                              elems[2],
                              elems[3]
                          )
                );
            }
        }
        for (int r = matrix.get_n(); r < this->rows; r++)
        {
            for (int c = 0; c < this->cols - 1; c++)
                this->set(r, c, _mm256_setzero_si256());

            /* lazy.... */
            switch (r % VECTOR_N)
            {
            case 1:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              1, 0, 0, 0
                         )
                );
                break;
            case 2:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 1ull << 32, 0, 0
                         )
                );
                break;
            case 3:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 1, 0, 0
                         )
                );
                break;
            case 4:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 0, 1ull << 32, 0
                         )
                );
                break;
            case 5:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 0, 1, 0
                         )
                );
                break;
            case 6:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 0, 0, 1ull << 32
                         )
                );
                break;
            case 7:
                this->set(r, this->cols - 1, _mm256_set_epi64x(
                              0, 0, 0, 1
                         )
                );
                break;
            }

        }
    }

    void mul_gamma(int r1, int r2, const GF_element &gamma)
    {
        /* here we do r1 first left to right and save the auxiliary gamma vectors.
         * then we permute the gamma vectors as required (note this->nmod here),
         * and do r2 left to right */
        /* do gamma multiplication during initialization?? */
        __m256i pac_gamma = _mm256_set_epi64x(
            gamma.get_repr(),
            gamma.get_repr(),
            gamma.get_repr(),
            gamma.get_repr()
        );
        // pac_gamma = [gamma^VECTOR_N]
        pac_gamma = global::F.wide_mul(pac_gamma, pac_gamma);
        pac_gamma = global::F.wide_mul(pac_gamma, pac_gamma);
        /* ugly */
        __m256i prod = _mm256_set_epi64x(
            1,
            gamma.get_repr(),
            (gamma*gamma).get_repr(),
            (gamma*gamma*gamma).get_repr()
        );
        std::vector<__m256i> coeffs(this->cols);
        /* first do left to right r1 */
        for (int col = 0; col < this->cols; col++)
        {
            /* already save them in reverse order here and permute,
             * values in reverse order, too*/
            coeffs[this->cols - 1 - col] = _mm256_permute4x64_epi64(
                prod,
                0x1B
            );
            __m256i elem = this->get(r1, col);
            elem = global::F.wide_mul(elem, prod);
            this->set(r1, col, elem);

            prod = global::F.wide_mul(prod, pac_gamma);
        }

        /* handle special permutations required in case not divisible by 4 */
        if (this->nmod)
        {
            __m256i idx;
            switch (this->nmod)
            {
            case 1:
                idx = _mm256_set_epi64x(0b000, 0b111, 0b110, 0b101);
                break;
            case 2:
                idx = _mm256_set_epi64x(0b001, 0b000, 0b111, 0b110);
                break;
            case 3:
                idx = _mm256_set_epi64x(0b010, 0b001, 0b000, 0b111);
                break;
            }

            for (int col = 0; col < this->cols - 1; col++)
                coeffs[col] = _mm256_permutex2var_epi64(
                    coeffs[col],
                    idx,
                    coeffs[col + 1]
                );

            switch (this->nmod)
            {
            case 1:
                coeffs[this->cols - 1] = _mm256_permute4x64_epi64(
                    coeffs[this->cols - 1],
                    0x00
                );
                break;
            case 2:
                coeffs[this->cols - 1] = _mm256_permute4x64_epi64(
                    coeffs[this->cols - 1],
                    0x40
                );
                break;
            case 3:
                coeffs[this->cols - 1] = _mm256_permute4x64_epi64(
                    coeffs[this->cols - 1],
                    0x90
                );
                break;
            }

        }

        /* and do r2 left to right */
        for (int col = 0; col < this->cols; col++)
        {
            this->set(r2, col,
                      global::F.wide_mul(
                          this->get(r2, col),
                          coeffs[col]
                      )
                );
        }
    }

    void swap_rows(int r1, int r2)
    {
        __m256i tmp;
        for (int c = 0; c < this->cols; c++)
        {
            tmp = this->get(r1, c);
            this->set(r1, c, this->get(r2, c));
            this->set(r2, c, tmp);
        }
    }

    void mul_row(int row, uint64_t v)
    {
        __m256i pack = _mm256_set_epi64x(
            v,
            v,
            v,
            v
        );
        for (int col = 0; col < this->cols; col++)
            this->set(row, col,
                      global::F.wide_mul(this->get(row, col), pack)
                );
    }

    /* subtract v times r1 from r2 */
    void row_op(int r1, int r2, uint64_t v)
    {
        __m256i pack = _mm256_set_epi64x(
            v,
            v,
            v,
            v
        );
        for (int col = 0; col < this->cols; col++)
        {
            __m256i tmp = global::F.wide_mul(this->get(r1, col), pack);

            this->set(r2, col,
                      _mm256_xor_si256(this->get(r2, col), tmp)
                );
        }
    }

    GF_element det()
    {
        uint64_t det = 0x1;
        for (int col = 0; col < this->cols; col++)
        {
            DET_LOOP(0);
            DET_LOOP(1);
            DET_LOOP(2);
            DET_LOOP(3);
        }
        return GF_element(det);
    }

    /* only used for testing */
    FMatrix unpack()
    {
        std::valarray<GF_element> unpacked(this->rows * this->rows);

        for (int row = 0; row < this->rows; row++)
        {
            for (int col = 0; col < this->cols; col++)
            {
                unpacked[row*this->rows + VECTOR_N*col + 0] =
                    GF_element(_mm256_extract_epi64(this->get(row, col), 3));
                unpacked[row*this->rows + VECTOR_N*col + 1] =
                    GF_element(_mm256_extract_epi64(this->get(row, col), 2));
                unpacked[row*this->rows + VECTOR_N*col + 2] =
                    GF_element(_mm256_extract_epi64(this->get(row, col), 1));
                unpacked[row*this->rows + VECTOR_N*col + 3] =
                    GF_element(_mm256_extract_epi64(this->get(row, col), 0));
            }
        }

        uint64_t n = this->rows;
        if (this->nmod)
            n -= VECTOR_N - this->nmod;
        return FMatrix(n, unpacked[std::gslice(0, {n,n}, {this->rows,1})]);
    }
};

#endif
