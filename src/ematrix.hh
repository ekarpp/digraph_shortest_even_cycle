/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef EMATRIX_H
#define EMATRIX_H

#include <valarray>

#include "extension.hh"
#include "matrix.hh"
#include "fmatrix.hh"

/* forward declare */
class FMatrix;

class EMatrix
{
private:
    Matrix<Extension_element> m;
    int n;

public:
    EMatrix(int n, std::valarray<Extension_element> m);
    EMatrix(Matrix<Extension_element> m);

    EMatrix operator+(const EMatrix &other) const;
    EMatrix operator-(const EMatrix &other) const;
    EMatrix operator*(const EMatrix &other) const;
    FMatrix project() const;

    const Matrix<Extension_element> &get_m() const { return this->m; }

    /* returns Per(this) - Det(this) as described in chapter 3
     * of the paper*/
    Extension_element per_m_det();

    Extension_element operator()(int row, int col) const
    {
        return this->m(row,col);
    }

    bool operator==(const EMatrix &other) const
    {
        return this->m == other.get_m();
    }

    bool operator!=(const EMatrix &other) const
    {
        return !(*this == other);
    }

    void set(int row, int col, Extension_element val)
    {
        this->m.set(row, col, val);
    }
};

#endif
