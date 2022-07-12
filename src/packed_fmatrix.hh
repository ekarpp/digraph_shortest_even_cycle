/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef P_FMATRIX_H
#define P_FMATRIX_H

#include <valarray>
#include <immintrin.h>

#include "gf.hh"
#include "global.hh"
#include "fmatrix.hh"

#define DET_LOOP(index)                                         \
    {                                                           \
        int r0 = 2*col + index;                                 \
        __m128i mx = this->get(r0, col);                        \
        int mxi = r0;                                           \
        uint64_t cmpmsk = 0b11 << (8*index);                    \
        for (int row = r0 + 1; row < this->rows; row++)         \
        {                                                       \
            uint64_t cmp = _mm_movemask_epi8(                   \
                _mm_cmpgt_epi16(                                \
                    this->get(row, col),                        \
                    mx                                          \
                    )                                           \
                );                                              \
            if (cmp & cmpmsk)                                   \
            {                                                   \
                mx = this->get(row, col);                       \
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
    std::valarray<__m128i> m;

    __m128i get(int row, int col)
    {
        return this->m[row*this->cols + col];
    }

    void set(int row, int col, __m128i v)
    {
        this->m[row*this->cols + col] = v;
    }

    __m128i gf_mul(__m128i a, __m128i b)
    {
        __m128i prod = _mm_unpacklo_epi64(
            _mm_clmulepi64_si128(a, b, 0x00),
            _mm_clmulepi64_si128(a, b, 0x11)
        );

        __m128i lo = _mm_shuffle_epi8(
            prod,
            _mm_set_epi32(
                0x80808080,
                0x80800908,
                0x80808080,
                0x80800100
            )
        );
        __m128i hi = _mm_shuffle_epi8(
            prod,
            _mm_set_epi32(
                0x80808080,
                0x80800B0A,
                0x80808080,
                0x80800302
            )
        );

        __m128i tmp = _mm_xor_si128(
            hi,
            _mm_xor_si128(
                _mm_srli_epi16(hi, 14),
                _mm_xor_si128(
                    _mm_srli_epi16(hi, 13),
                    _mm_srli_epi16(hi, 11)
                )
            )
        );

        __m128i rem = _mm_xor_si128(
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

    void mul_gamma(int r1, int r2, const GF_element &gamma)
    {
        __m128i pac_gamma = _mm_set_epi64x(
            gamma.get_repr(),
            gamma.get_repr()
        );
        pac_gamma = this->gf_mul(pac_gamma, pac_gamma);
        __m128i prod = _mm_set_epi64x(
            1,
            gamma.get_repr()
        );
        /* first do left to right r1 */
        for (int col = 0; col < this->cols; col++)
        {
            __m128i elem = this->get(r1, col);
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
            __m128i elem = this->get(r2, col);
            elem = this->gf_mul(elem, prod);
            this->set(r2, col, elem);

            prod = this->gf_mul(prod, pac_gamma);
        }
    }

    void swap_rows(int r1, int r2)
    {
        __m128i tmp;
        for (int c = 0; c < this->cols; c++)
        {
            tmp = this->get(r1, c);
            this->set(r1, c, this->get(r2, c));
            this->set(r2, c, tmp);
        }
    }

    void mul_row(int row, uint64_t v)
    {
        __m128i pack = _mm_set_epi64x(
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
        __m128i pack = _mm_set_epi64x(
            v,
            v
        );
        for (int col = 0; col < this->cols; col++)
        {
            __m128i tmp = this->gf_mul(this->get(r1, col), pack);

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
};

#endif
