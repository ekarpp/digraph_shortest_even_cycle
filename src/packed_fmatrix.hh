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
        uint64_t v = this->m[row*this->cols + col];
        /* can skip right shift?? */
        return (col%2 == 0)
            ? v >> 32
            : v & 0xFFFF;
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
            elem = global::F.packed_clmul(elem, prod);
            elem = global::F.packed_rem(elem);
            this->m[r1*this->cols + col] = elem;

            prod = global::F.packed_clmul(prod, pac_gamma);
            prod = global::F.packed_rem(prod);
        }
        /* swap prod around */
        prod = (prod << 32) | (prod >> 32);
        /* pack inverse of gamma */
        pac_gamma = global::F.ext_euclid(gamma.get_repr());
        pac_gamma |= pac_gamma << 32;
        /* fix edge case of odd n */
        if (this->rows % 2)
        {
            prod = global::F.packed_clmul(prod, pac_gamma);
            prod = global::F.packed_rem(prod);
        }
        /* now do left to right r2 */
        for (int col = 0; col < this->cols; col++)
        {
            uint64_t elem = this->m[r2*this->cols + col];
            elem = global::F.packed_clmul(elem, prod);
            elem = global::F.packed_rem(elem);
            this->m[r2*this->cols + col] = elem;

            prod = global::F.packed_clmul(prod, pac_gamma);
            prod = global::F.packed_rem(prod);
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
            this->m[idx] = global::F.packed_clmul(this->m[idx], v);
            this->m[idx] = global::F.packed_rem(this->m[idx]);
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
            uint64_t tmp = global::F.packed_clmul(this->m[idx1], v);
            tmp = global::F.packed_rem(tmp);

            this->m[idx2] ^= tmp;
        }
    }

    GF_element det()
    {
        uint64_t det = 0x1;
        for (int col = 0; col < this->cols; col++)
        {
            for (int i = 0; i < 2; i++)
            {
                int c0 = 2*col + i;
                uint64_t mx = this->get(c0, c0);
                int mxi = c0;

                for (int row = c0 + 1; row < this->rows; row++)
                {
                    if (this->get(row, c0) > mx)
                    {
                        mx = this->get(row, c0);
                        mxi = row;
                    }
                }

                if (mx == 0)
                    return global::F.zero();

                if (mxi != c0)
                    this->swap_rows(mxi, c0);

                det = global::F.clmul(det, mx);
                det = global::F.rem(det);
                mx = global::F.ext_euclid(mx);

                this->mul_row(c0, mx);
                for (int row = c0 + 1; row < this->rows; row++)
                    this->row_op(c0, row, this->get(row, c0));
            }
        }
        return GF_element(det);
    }
};

#endif
