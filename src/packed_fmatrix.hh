/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef P_FMATRIX_H
#define P_FMATRIX_H

#include <valarray>
#include <immintrin.h>

#include "gf.hh"
#include "global.hh"
#include "fmatrix.hh"

typedef long long int long2_t __attribute__ ((vector_size (16)));

#define DET_LOOP(index)                                         \
    {                                                           \
        int r0 = 2*col + index;                                 \
        long2_t mx = this->get(r0, col);                        \
        int mxi = r0;                                           \
        uint64_t cmpmsk = 0b11 << (8*(1-index));                \
        for (int row = r0 + 1; row < this->rows; row++)         \
        {                                                       \
            uint64_t cmp = _mm_movemask_epi8(                   \
                _mm_cmpgt_epi32(                                \
                    this->get(row, col),                        \
                    mx                                          \
                    )                                           \
                );                                              \
            if (cmp & cmpmsk)                                   \
            {                                                   \
                mx = this->get(row,col);                        \
                mxi = row;                                      \
            }                                                   \
        }                                                       \
        uint64_t mx_ext = _mm_extract_epi64(mx, 1 - index);     \
        if (mx_ext == 0)                                        \
            return global::F.zero();                            \
        if (mxi != r0)                                          \
            this->swap_rows(mxi, r0);                           \
        det = global::F.rem(                                    \
            global::F.clmul(det, mx_ext)                        \
            );                                                  \
        mx_ext = global::F.ext_euclid(mx_ext);                  \
        this->mul_row(r0, mx_ext);                              \
        for (int row = r0 + 1; row < this->rows; row++)         \
        {                                                       \
            uint64_t val = _mm_extract_epi64(                   \
                this->get(row, col),                            \
                1 - index                                       \
                );                                              \
            this->row_op(r0, row, val);                         \
        }                                                       \
    }

class Packed_FMatrix
{
private:
    int rows;
    int cols;
    std::valarray<long2_t> m;

    long2_t get(int row, int col)
    {
        return this->m[row*this->cols + col];
    }

    void set(int row, int col, long2_t v)
    {
        this->m[row*this->cols + col] = v;
    }

    long2_t gf_mul(long2_t a, long2_t b)
    {
        long2_t prod = _mm_unpacklo_epi64(
            _mm_clmulepi64_si128(a, b, 0x00),
            _mm_clmulepi64_si128(a, b, 0x11)
        );

        long2_t lo = _mm_shuffle_epi8(
            prod,
            _mm_set_epi32(
                0x80808080,
                0x80800908,
                0x80808080,
                0x80800100
            )
        );
        long2_t hi = _mm_shuffle_epi8(
            prod,
            _mm_set_epi32(
                0x80808080,
                0x80800B0A,
                0x80808080,
                0x80800302
            )
        );

        long2_t tmp = _mm_xor_si128(
            hi,
            _mm_xor_si128(
                _mm_srli_epi16(hi, 14),
                _mm_xor_si128(
                    _mm_srli_epi16(hi, 13),
                    _mm_srli_epi16(hi, 11)
                )
            )
        );

        long2_t rem = _mm_xor_si128(
            tmp,
            _mm_xor_si128(
                _mm_slli_epi16(tmp, 2),
                _mm_xor_si128(
                    _mm_slli_epi16(tmp, 3),
                    _mm_slli_epi16(tmp, 5)
                )
            )
        );

        return _mm_xor_si128(rem, lo);
    }

public:
    Packed_FMatrix(
        const FMatrix &matrix
    )
    {
        this->rows = matrix.get_n();
        this->cols = this->rows / 2 + (this->rows % 2);

        this->m.resize(this->rows * this->cols);

        for (int r = 0; r < this->rows; r++)
        {
            for (int c = 0; c < this->rows / 2; c++)
                this->set(r, c,
                          _mm_set_epi64x(
                              matrix(r, 2*c).get_repr(),
                              matrix(r, 2*c + 1).get_repr()
                          )
                    );


            if (this->rows % 2)
                this->set(r, this->cols - 1,
                          _mm_set_epi64x(
                              matrix(r, this->rows - 1).get_repr(),
                              0
                          )
                    );
        }
    }

    void handle_gamma_zero(int r1, int r2)
    {
        long2_t mask = _mm_set_epi32(
            0xFFFFFFFFull,
            0xFFFFFFFFull,
            0x0,
            0x0
        );
        this->set(r1, 0,
                  _mm_and_si128(
                      this->get(r1, 0),
                      mask
                  )
            );
        if (this->rows % 2 == 0)
            mask = _mm_shuffle_epi32(
                mask,
                0b00011011
            );

        this->set(r2, this->cols - 1,
                  _mm_and_si128(
                      this->get(r2, this->cols - 1),
                      mask
                   )
            );

        for (int col = 1; col < this->cols - 1; col++)
        {
            this->set(r1, col, _mm_setzero_si128());
            this->set(r2, col, _mm_setzero_si128());
        }
        this->set(r2, 0, _mm_setzero_si128());
        this->set(r1, this->cols - 1, _mm_setzero_si128());
    }

    void mul_gamma(int r1, int r2, const GF_element &gamma)
    {
        /* alternative approach: do multiplications only in one direction
         * and then in the other direction use shuffle and shift them around
         * as neccessary (note the uneven row case).
         * this should be fewer multiplications (and instructions)
         * and don't need to separately handle the case of gamma being zero. */
        /* do gamma multiplication during initialization?? */
        long2_t pac_gamma = _mm_set_epi64x(
            gamma.get_repr(),
            gamma.get_repr()
        );
        pac_gamma = this->gf_mul(pac_gamma, pac_gamma);
        long2_t prod = _mm_set_epi64x(
            1,
            gamma.get_repr()
        );
        /* first do left to right r1 */
        for (int col = 0; col < this->cols; col++)
        {
            long2_t elem = this->get(r1, col);
            elem = this->gf_mul(elem, prod);
            this->set(r1, col, elem);

            prod = this->gf_mul(prod, pac_gamma);
        }
        /* swap prod around */
        prod = _mm_shuffle_epi32(prod, 0b01001110);

        /* pack inverse of gamma */
        uint64_t gamma_inv = gamma.inv().get_repr();
        pac_gamma = _mm_set_epi64x(
            gamma_inv,
            gamma_inv
        );
        /* fix edge case of odd n, degree is one too high */
        if (this->rows % 2)
            prod = this->gf_mul(prod, pac_gamma);
        pac_gamma = this->gf_mul(pac_gamma, pac_gamma);
        /* final iteration did one extra multiply, revert it */
        prod = this->gf_mul(prod, pac_gamma);
        /* now do left to right r2 */
        for (int col = 0; col < this->cols; col++)
        {
            long2_t elem = this->get(r2, col);
            elem = this->gf_mul(elem, prod);
            this->set(r2, col, elem);

            prod = this->gf_mul(prod, pac_gamma);
        }
    }

    void swap_rows(int r1, int r2)
    {
        long2_t tmp;
        for (int c = 0; c < this->cols; c++)
        {
            tmp = this->get(r1, c);
            this->set(r1, c, this->get(r2, c));
            this->set(r2, c, tmp);
        }
    }

    void mul_row(int row, uint64_t v)
    {
        long2_t pack = _mm_set_epi64x(
            v,
            v
        );
        for (int col = 0; col < this->cols; col++)
            this->set(row, col,
                      this->gf_mul(this->get(row, col), pack)
                );
    }

    /* subtract v times r1 from r2 */
    void row_op(int r1, int r2, uint64_t v)
    {
        long2_t pack = _mm_set_epi64x(
            v,
            v
        );
        for (int col = 0; col < this->cols; col++)
        {
            long2_t tmp = this->gf_mul(this->get(r1, col), pack);

            this->set(r2, col,
                      _mm_xor_si128(this->get(r2, col), tmp)
                );
        }
    }

    GF_element det()
    {
        uint64_t det = 0x1;
        for (int col = 0; col < this->cols - (this->rows%2); col++)
        {
            DET_LOOP(0);
            DET_LOOP(1);
        }
        if (this->rows % 2)
        {
            uint64_t mx_ext = _mm_extract_epi64(
                this->get(this->rows - 1, this->cols - 1),
                1
            );
            det = global::F.rem(
                global::F.clmul(det, mx_ext)
            );
        }
        return GF_element(det);
    }

    /* only used for testing */
    FMatrix unpack()
    {
        std::valarray<GF_element> unpacked(this->rows * this->rows);

        for (int row = 0; row < this->rows; row++)
        {
            for (int col = 0; col < this->cols - (this->rows%2); col++)
            {
                unpacked[row*this->rows + 2*col] =
                    GF_element(_mm_extract_epi64(this->get(row, col), 1));
                unpacked[row*this->rows + 2*col + 1] =
                    GF_element(_mm_extract_epi64(this->get(row, col), 0));
            }
            if (this->rows%2)
            {
                int col = this->cols - 1;
                unpacked[row*this->rows + 2*col] =
                    GF_element(_mm_extract_epi64(this->get(row, col), 1));
            }
        }
        return FMatrix(this->rows, unpacked);
    }
};

#endif
