#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

/* parent class for square matrices with different elements
   contains the usual arithmetic operations */

template <typename E>
class Matrix
{
    using _et = E;

private:
    std::vector<std::vector<_et>> m;
    int n;
    _et one;

public:
    Matrix(int n, _et one, std::vector<std::vector<_et>> &m)
    {
        this->n = n;
        this->one = one;
        this->m(n, std::vector<_et>(n));
    }

    Matrix &operator+(const Matrix &other)
    {
        std::vector<std::vector<_et>> sum(n, std::vector<_et>(n));

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; i++)
                sum[i][j] = *this[i][j] + other[i][j];

        return Matrix(this->n, this->one, sum);
    }

    Matrix &operator-(const Matrix &other)
    {
        std::vector<std::vector<_et>> sum(n, std::vector<_et>(n));

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; i++)
                sum[i][j] = *this[i][j] - other[i][j];

        return Matrix(this->n, this->one, sum);
    }

    Matrix &operator*(const Matrix &other)
    {
        std::vector<std::vector<_et>> prod(n, std::vector<_et>(n));

        for (int i = 0; i < this->n; i++)
        {
            for (int j = 0; j < this->n; j++)
            {
                prod[i][j] = this->one;
                for (int k = 0; k < this->n; k++)
                    prod[i][j] = prod[i][j] + *this[i][j] * other[i][j];
            }
        }

        return Matrix(this->n, this->one, prod);
    }

    std::vector<_et> &operator[](int i)
    {
        /* todo: error check, append zeros? */
        return this->m[i];
    }

    int get_n() { return this->n; }
};

#endif
