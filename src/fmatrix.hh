#ifndef FMATRIX_H
#define FMATRIX_H

#include <vector>

#include "gf.hh"
#include "matrix.hh"

class FMatrix
{
private:
    Matrix<GF_element> m;
    int n;

public:
    FMatrix(int n, std::vector<std::vector<GF_element>> m);
    FMatrix(Matrix<GF_element> m);

    FMatrix operator+(const FMatrix &other) const;

    const Matrix<GF_element> &get_m() const { return this->m; }

    std::vector<GF_element> operator[](int i)
    {
        return this->m[i];
    }
};

#endif
