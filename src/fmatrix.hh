#ifndef FMATRIX_H
#define FMATRIX_H

#include <vector>

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
    FMatrix(int n, std::vector<std::vector<GF_element>> m);
    FMatrix(Matrix<GF_element> m);

    FMatrix operator+(const FMatrix &other) const;
    FMatrix operator-(const FMatrix &other) const;
    FMatrix operator*(const FMatrix &other) const;

    EMatrix lift() const;

    const Matrix<GF_element> &get_m() const { return this->m; }

    std::vector<GF_element> operator[](int i) const
    {
        return this->m[i];
    }

    GF_element operator()(int row, int col) const
    {
        /* append identity matrix in case out of bounds.
         * potentially dangerous */
        if (row >= this->n || col >= this->n)
            return (row == col) ? global::F.one() : global::F.zero();

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
};

#endif
