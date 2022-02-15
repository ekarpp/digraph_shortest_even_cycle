#include <valarray>

#include "global.hh"
#include "fmatrix.hh"
#include "ematrix.hh"
#include "extension.hh"

using namespace std;

FMatrix::FMatrix(int n, valarray<GF_element> m): m(n, global::F.one(), m)
{
    this->n = n;
}

FMatrix::FMatrix(Matrix<GF_element> m)
{
    this->n = m.get_n();
    this->m = m;
}

FMatrix FMatrix::operator+(const FMatrix &other) const
{
    return FMatrix(this->m + other.get_m());
}

FMatrix FMatrix::operator-(const FMatrix &other) const
{
    return FMatrix(this->m - other.get_m());
}

FMatrix FMatrix::operator*(const FMatrix &other) const
{
    return FMatrix(this->m * other.get_m());
}

EMatrix FMatrix::lift() const
{
    valarray<Extension_element> lifted(this->n * this->n);

    for (int x = 0; x < this->n; x++)
    {
        for (int y = 0; y < this->n; y++)
            lifted[x*this->n + y] = this->operator()(x,y).lift();
    }

    return EMatrix(this->n, lifted);
}

/* PA = LU, L lower triangular with ones on diagonal
 * and U upper triangular. P permutation matrix as valarray.
   modifies the object it is called on to L and U in single matrix. */
valarray<int> FMatrix::lup()
{
    valarray<int> P(1,1);

    return P;
}

GF_element FMatrix::det() const
{
    FMatrix LU = this->copy();
    valarray<int> P = LU.lup();
    GF_element det = global::F.one();
    /* PA = LU => det(A) = det(P^T)det(L)det(U)
     * in current implementation L has diagonal of one
     * thus det(L) = 1. det(U) is the product of its diagonal
     * elements. det(P)=det(P^T) is the sgn of P. */
    for (int i = 0; i < this->n; i++)
        det *= LU(i,i);
    /* check sgn of P. write func in util. its used elsewhere also. */
    return det;
}

FMatrix FMatrix::copy() const
{
    valarray<GF_element> m(this->n * this->n);

    for (int row = 0; row < this->n; row++)
        for (int col = 0; col < this->n; col++)
            m[row*n + col] = this->m(row,col);

    return FMatrix(this->n, m);
}
