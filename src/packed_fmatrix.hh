/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef P_FMATRIX_H
#define P_FMATRIX_H

#include <valarray>

#include "gf.hh"
#include "global.hh"
#include "fmatrix.hh"

class Packed_FMatrix
{
private:
    int rows;
    int cols;
    std::valarray<uint64_t> m;

    uint64_t get(int row, int col)
    {
        uint64_t v = this->m[row*this->cols + col / 2];
        /* can skip right shift?? */
        return (col%2 == 0)
            ? v >> 32
            : v & 0xFFFF;
    }

    uint64_t gf_mul(uint64_t a, uint64_t b)
    {
        return global::F.packed_rem(
            global::F.packed_clmul(a, b)
        );
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
            {
                this->m[r*this->cols + c] = matrix(r, 2*c).get_repr() << 32;
                this->m[r*this->cols + c] |= matrix(r, 2*c + 1).get_repr();
            }
            if (this->rows % 2)
                this->m[r*this->cols + this->rows/2] =
                    matrix(r, this->rows - 1).get_repr() << 32;
        }
    }

    void mul_gamma(int r1, int r2, const GF_element &gamma)
    {
        uint64_t pac_gamma = gamma.get_repr() | (gamma.get_repr() << 32);
        uint64_t prod = gamma.get_repr() | (1ull << 32);
        /* first do left to right r1 */
        for (int col = 0; col < this->cols; col++)
        {
            uint64_t elem = this->m[r1*this->cols + col];
            elem = this->gf_mul(elem, prod);
            this->m[r1*this->cols + col] = elem;

            prod = this->gf_mul(prod, pac_gamma);
        }
        /* swap prod around */
        prod = (prod << 32) | (prod >> 32);
        /* pack inverse of gamma */
        pac_gamma = global::F.ext_euclid(gamma.get_repr());
        pac_gamma |= pac_gamma << 32;
        /* fix edge case of odd n */
        if (this->rows % 2)
            prod = this->gf_mul(prod, pac_gamma);

        /* now do left to right r2 */
        for (int col = 0; col < this->cols; col++)
        {
            uint64_t elem = this->m[r2*this->cols + col];
            elem = this->gf_mul(elem, prod);
            this->m[r2*this->cols + col] = elem;

            prod = this->gf_mul(prod, pac_gamma);
        }
    }

    void swap_rows(int r1, int r2)
    {
        uint64_t tmp;
        for (int c = 0; c < this->cols; c++)
        {
            int idx1 = r1*this->cols + c;
            int idx2 = r2*this->cols + c;
            tmp = this->m[idx1];
            this->m[idx1] = this->m[idx2];
            this->m[idx2] = tmp;
        }
    }

    void mul_row(int row, uint64_t v)
    {
        v |= v << 32;
        for (int col = 0; col < this->cols; col++)
        {
            int idx = row*this->cols + col;
            this->m[idx] = this->gf_mul(this->m[idx], v);
        }
    }

    /* subtract v times r1 from r2 */
    void row_op(int r1, int r2, uint64_t v)
    {
        /* pack v */
        v |= v << 32;
        for (int col = 0; col < this->cols; col++)
        {
            int idx1 = r1*this->cols + col;
            int idx2 = r2*this->cols + col;
            uint64_t tmp = this->gf_mul(this->m[idx1], v);

            this->m[idx2] ^= tmp;
        }
    }

    GF_element det()
    {
        uint64_t det = 0x1;
        for (int col = 0; col < this->rows; col++)
        {
            //int c0 = 2*col + i;
            uint64_t mx = this->get(col, col);
            int mxi = col;

            for (int row = col + 1; row < this->rows; row++)
            {
                if (this->get(row, col) > mx)
                {
                    mx = this->get(row, col);
                    mxi = row;
                }
            }

            if (mx == 0)
                return global::F.zero();

            if (mxi != col)
                this->swap_rows(mxi, col);

            det = global::F.rem(
                global::F.clmul(det, mx)
            );
            mx = global::F.ext_euclid(mx);

            this->mul_row(col, mx);
            for (int row = col + 1; row < this->rows; row++)
                this->row_op(col, row, this->get(row, col));
        }
        return GF_element(det);
    }
};

#endif
