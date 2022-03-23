/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef FMATRIX_H
#define FMATRIX_H

#include <iostream>
#include <bitset>
#include <vector>
#include <valarray>

#include "gf.hh"
#include "ematrix.hh"
#include "matrix.hh"
#include "polynomial.hh"

/* forward declare */
class EMatrix;

class FMatrix
{
private:
    Matrix<GF_element> m;
    int n;

public:
    /* for graph.cc */
    FMatrix() {};
    FMatrix(int n, std::valarray<GF_element> m);
    FMatrix(Matrix<GF_element> m);

    FMatrix operator+(const FMatrix &other) const;
    FMatrix operator-(const FMatrix &other) const;
    FMatrix operator*(const FMatrix &other) const;

    /* return pcc_{n-1} of the matrix we get when we
     * multiply the diagonal of this matrix by e */
    GF_element pcc(GF_element e) const;

    EMatrix lift() const;

    /* multiply diagonal by e. merge this with lift, so
     * that only one new copy is created? lift gets always
     * called after this */
    FMatrix mul_diag(GF_element e) const;

    /* returns a permutation matrix P and modifies the object called on
     * such that P*obj = obj.lup(). where obj.lup() contains L on the lower
     * triangle and U on the upper triangle such that P*obj = L*U. note that
     * L has ones on diagonal, thus L and U can be stored in one matrix. */
    std::vector<int> lup(int depth);

    /* uses lup(). modifies the object it is called on. */
    GF_element det();

    /* det of the matrix we get when r1 is multiplied by monomials
     * (1,r,..,r^(n-1)) and r2 by monomials (r^(n-1),..,r,1) */
    Polynomial pdet(int r1, int r2) const;

    /* returns a copy of this */
    FMatrix copy() const;

    const Matrix<GF_element> &get_m() const { return this->m; }

    void mul(int row, int col, GF_element v)
    {
        this->m.mul(row, col, v);
    }

    void mul_row(int row, GF_element v)
    {
        this->m.mul_row(row, v);
    }

    /* subtract v times r1 from r2 */
    void row_op(int r1, int r2, GF_element v)
    {
        this->m.row_op(r1, r2, v);
    }

    GF_element operator()(int row, int col) const
    {
        return this->m(row,col);
    }

    bool operator==(const FMatrix &other) const
    {
        return this->m == other.get_m();
    }

    bool operator!=(const FMatrix &other) const
    {
        return !(*this == other);
    }

    void set(int row, int col, GF_element val)
    {
        this->m.set(row, col, val);
    }

    /* swap rows r1 and r2 starting from column idx */
    void swap_rows(int r1, int r2, int idx = 0)
    {
        GF_element tmp;
        for (int col = idx; col < this->n; col++)
        {
            tmp = this->m(r1, col);
            this->set(r1, col, this->m(r2,col));
            this->set(r2, col, tmp);
        }
    }

    std::valarray<GF_element> slice(int start, int size, int stride) const
    {
        return this->m.slice(start, size, stride);
    }

    void print() const
    {
        for (int row = 0; row < this->n; row++)
        {
            for (int col = 0; col < this->n; col++)
                std::cout << std::bitset<8>(this->m(row, col).get_repr()) << " ";
            std::cout << std::endl;
        }
    }
};

#endif
