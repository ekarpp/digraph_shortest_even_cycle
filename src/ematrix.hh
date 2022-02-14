#ifndef EMATRIX_H
#define EMATRIX_H

#include <vector>

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
    EMatrix(int n, std::vector<std::vector<Extension_element>> m);
    EMatrix(Matrix<Extension_element> m);

    EMatrix operator+(const EMatrix &other) const;
    EMatrix operator-(const EMatrix &other) const;
    EMatrix operator*(const EMatrix &other) const;
    FMatrix project() const;

    const Matrix<Extension_element> &get_m() const { return this->m; }

    std::vector<Extension_element> operator[](int i)
    {
        return this->m[i];
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
