/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <valarray>
#include <list>

#include "polynomial.hh"
#include "global.hh"
#include "ematrix.hh"
#include "fmatrix.hh"
#include "gf.hh"

using namespace std;

EMatrix::EMatrix(int n, valarray<Extension_element> m): m(n, global::E.one(), m)
{
    this->n = n;
}

EMatrix::EMatrix(Matrix<Extension_element> m)
{
    this->n = m.get_n();
    this->m = m;
}

EMatrix EMatrix::operator+(const EMatrix &other) const
{
    return EMatrix(this->m + other.get_m());
}

EMatrix EMatrix::operator-(const EMatrix &other) const
{
    return EMatrix(this->m - other.get_m());
}

EMatrix EMatrix::operator*(const EMatrix &other) const
{
    return EMatrix(this->m * other.get_m());
}

EMatrix EMatrix::copy() const
{
    valarray<Extension_element> m(this->n * this->n);

    for (int row = 0; row < this->n; row++)
        for (int col = 0; col < this->n; col++)
            m[row*this->n + col] = this->m(row, col);

    return EMatrix(this->n, m);
}

FMatrix EMatrix::project() const
{
    valarray<GF_element> proj(this->n * this->n);

    for (int x = 0; x < this->n; x++)
    {
        for (int y = 0; y < this->n; y++)
            proj[x*this->n + y] = this->operator()(x,y).project();
    }

    return FMatrix(this->n, proj);
}

/* returns Per(this) - Det(this) as described in chapter 3
 * of the paper*/
Extension_element EMatrix::per_m_det()
{
    Extension_element acc = global::E.zero();

    /* marked rows */
    valarray<bool> rows(false, this->n);
    /* odd elements at (odd[i], i). if odd[i] = -1 then
     * column i has only even elements at unmarked rows */
    valarray<int> odd(this->n);
    /* columns that have only even elements at unmarked rows */
    list<int> cols;

    for (int j = 0; j < this->n; j++)
    {
        int i1;
        for (i1 = 0; j < this->n; j++)
        {
            if (rows[i1])
                continue;
            /* transpose? */
            if (!this->operator()(i1, j).is_even())
            {
                acc += this->row_op(i1, j);
                rows[i1] = true;
                odd[j] = i1;
                break;
            }
        }
        if (i1 == this->n)
        {
            odd[j] = -1;
            cols.push_front(j);
        }
    }

    Extension_element det = global::E.zero();
    /* if more than two unmarked columns, det and per
     * of the final matrix is zero because in characteristic
     * 2 even*even = 0 */
    if (cols.size() <= 1)
    {
        if (cols.size() == 1)
        {
            int row;
            /* find unmarked row */
            for (row = 0; row < this->n; row++)
                if (!rows[row])
                    break;
            /* unmarked column */
            int col = cols.front();
            odd[col] = row;
        }

        /* permanent is the product of the odd and maybe
         * one even element at the crossing of unmarked row
         * and column */
        int swaps = 0;
        valarray<bool> swapped(false, this->n);
        Extension_element per = global::E.one();
        for (int col = 0; col < this->n; col++)
        {
            int row = odd[col];
            per *= this->operator()(row, col);
            /* works? */
            if (row != col)
                swaps++;
        }
        acc += per;

        /* permutation sign */
        /* can just skip this and not just add per to acc */
        if (swaps % 2 == 1)
            /* unary - ? */
            det = global::E.zero() - per;
        else
            det = per;
    }

    return acc - det;
}

/* make all elements in row j even except for (i1,j)
 * return accumulator */
Extension_element EMatrix::row_op(int i1, int j)
{
    Extension_element acc = global::E.zero();
    const Extension_element sigma = this->operator()(i1, j);
    for (int i2 = 0; i2 < this->n; i2++)
    {
        if (i2 == i1)
            continue;
        if (!this->operator()(i2, j).is_even())
        {
            const Extension_element v = this->operator()(i2, j);
            const Extension_element t = util::tau(sigma, v);

            /* M'' in the paper. Modify this as M' */
            EMatrix mpp = this->copy();

            for (int col = 0; col < this->n; col++)
            {
                this->set(i2, col,
                          this->operator()(i2, col) - t*this->operator()(i1, col)
                );
                mpp.set(i2, col, t * this->operator()(i1, col));
            }
            /* project creates new copy here */
            const Polynomial p = mpp.project().pdet(i1, i2);
            GF_element sum = global::F.zero();
            for (int i = 0; i < this->n - 1; i++)
                sum += p[i];
            Extension_element per = sum.lift() + sum.lift();
            acc += per;
        }
    }
    return acc;
}
