#ifndef FMATRIX_H
#define FMATRIX_H

#include <vector>
#include <valarray>

#include "gf.hh"
#include "ematrix.hh"
#include "matrix.hh"

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

    EMatrix lift() const;

    /* returns a permutation matrix P and modifies the object called on
     * such that P*obj = obj.lup(). where obj.lup() contains L on the lower
     * triangle and U on the upper triangle such that P*obj = L*U. note that
     * L has ones on diagonal, thus L and U can be stored in one matrix. */
    std::vector<int> lup(int depth);

    /* uses lup(). calls it on copy of the object it is on called on */
    GF_element det() const;

    /* returns a copy of this */
    FMatrix copy() const;

    const Matrix<GF_element> &get_m() const { return this->m; }

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

};

#endif
