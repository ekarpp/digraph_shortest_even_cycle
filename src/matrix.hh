#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

/* parent class for square matrices with different elements
   contains the usual arithmetic operations */

template <typename E>
class Matrix
{
private:
    std::vector<std::vector<E>> m;
    int n;
    E one;

public:
    Matrix() {}

    Matrix(int n, E one, std::vector<std::vector<E>> &matrix)
    {
        this->n = n;
        this->one = one;
        std::vector<std::vector<E>> v(n, std::vector<E>(n));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                v[i][j] = matrix[i][j];
        this->m = v;
    }

    Matrix &operator+=(const Matrix &other)
    {
        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                this[i][j] = this[i][j] + other[i][j];

        return *this;
    }

    Matrix operator+(const Matrix &other) const
    {
        std::vector<std::vector<E>> sum(n, std::vector<E>(n));

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                sum[i][j] = (*this)[i][j] + other[i][j];

        return Matrix(this->n, this->one, sum);
    }

    Matrix operator-(const Matrix &other)
    {
        std::vector<std::vector<E>> sum(n, std::vector<E>(n));

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                sum[i][j] = this[i][j] - other[i][j];

        return Matrix(this->n, this->one, sum);
    }

    Matrix operator*(const Matrix &other)
    {
        std::vector<std::vector<E>> prod(n, std::vector<E>(n));

        for (int i = 0; i < this->n; i++)
        {
            for (int j = 0; j < this->n; j++)
            {
                prod[i][j] = this->one;
                for (int k = 0; k < this->n; k++)
                    prod[i][j] = prod[i][j] + this[i][j] * other[i][j];
            }
        }

        return Matrix(this->n, this->one, prod);
    }

    std::vector<E> operator[](int i)
    {
        /* todo: error check, append zeros? */
        return this->m[i];
    }

    int get_n() { return this->n; }
};

#endif
