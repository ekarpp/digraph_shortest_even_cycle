/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef MATRIX_H
#define MATRIX_H

#include <valarray>

/* parent class for square matrices with different elements
   contains the usual arithmetic operations */

template <typename E>
class Matrix
{
private:
    std::valarray<E> m;
    int n;
    E one;
    E zero;

public:
    Matrix() {}

    Matrix(int n, E one, std::valarray<E> &matrix)
    {
        this->n = n;
        this->one = one;
        this->zero = one - one;
        std::valarray<E> v(n * n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                v[i*n + j] = matrix[i*n + j];
        this->m = v;
    }

    Matrix &operator+=(const Matrix &other)
    {
        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                this->set(i, j, this->operator()(i,j) + other(i,j));

        return *this;
    }

    Matrix operator+(const Matrix &other) const
    {
        std::valarray<E> sum(this->n * this->n);

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                sum[i*this->n + j] = this->operator()(i,j) + other(i,j);

        return Matrix(this->n, this->one, sum);
    }

    Matrix operator-(const Matrix &other) const
    {
        std::valarray<E> sum(this->n * this->n);

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                sum[i*this->n + j] = this->operator()(i,j) - other(i,j);

        return Matrix(this->n, this->one, sum);
    }

    Matrix operator*(const Matrix &other) const
    {
        std::valarray<E> prod(this->n * this->n);

        for (int i = 0; i < this->n; i++)
        {
            for (int j = 0; j < this->n; j++)
            {
                prod[i*this->n + j] = this->zero;
                for (int k = 0; k < this->n; k++)
                    prod[i*this->n + j] += this->operator()(i,k) * other(k,j);
            }
        }

        return Matrix(this->n, this->one, prod);
    }

    const E &operator()(int row, int col) const
    {
        /* append identity matrix in case out of bounds.
         * potentially dangerous */
        /* commented for now, no use */
        //if (row >= this->n || col >= this->n)
        //    return (row == col) ? this->one : this->zero;
        return this->m[row*n + col];
    }

    int get_n() const { return this->n; }

    bool operator==(const Matrix &other) const
    {
        if (this->n != other.get_n())
            return false;

        for (int i = 0; i < this->n; i++)
            for (int j = 0; j < this->n; j++)
                if (this->operator()(i,j) != other(i,j))
                    return false;

        return true;
    }

    bool operator!=(const Matrix &other) const
    {
        return !(*this == other);
    }

    /* set m[row][col] = val. use this to avoid const issues with [] op */
    void set(int row, int col, E val)
    {
        this->m[row*n + col] = val;
    }

    /* multiply element at row,col by val */
    void mul(int row, int col, const E &val)
    {
        this->m[row*n + col] *= val;
    }

    /* multiply row by val */
    void mul_row(int row, const E &val)
    {
        for (int col = 0; col < this->n; col++)
            this->m[row*n + col] *= val;
    }

    /* subtract v times r1 from r2 */
    void row_op(int r1, int r2, const E &v)
    {
        for (int col = 0; col < this->n; col++)
            this->m[r2*n + col] -= v*this->operator()(r1,col);
    }

    /* subtract val from element at row,col */
    void sub(int row, int col, const E &val)
    {
        this->m[row*n + col] -= val;
    }

    std::valarray<E> slice(int start, int size, int stride) const
    {
        return this->m[std::slice(start, size, stride)];
    }
};

#endif
