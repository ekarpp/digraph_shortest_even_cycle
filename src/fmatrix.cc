#include <vector>

#include "global.hh"
#include "fmatrix.hh"
#include "ematrix.hh"
#include "extension.hh"

using namespace std;

FMatrix::FMatrix(int n, vector<vector<GF_element>> m): m(n, global::F.one(), m)
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
    vector<vector<Extension_element>> lifted(this->n, vector<Extension_element>(this->n));

    for (int x = 0; x < this->n; x++)
    {
        const vector<GF_element> row = this->m[x];
        for (int y = 0; y < this->n; y++)
            lifted[x][y] = row[y].lift();
    }

    return EMatrix(this->n, lifted);
}

/* PA = LU, L lower triangular with ones on diagonal
 * and U upper triangular. P permutation matrix as vector.
   modifies the object it is called on to L and U in single matrix. */
vector<int> FMatrix::lup()
{
    vector<int> P(1,1);

    return P;
}

GF_element FMatrix::det() const
{
    FMatrix LU = this->copy();
    vector<int> P = LU.lup();
    GF_element det = global::F.one();
    /* PA = LU => det(A) = det(P^T)det(L)det(U)
     * in current implementation L has diagonal of one
     * thus det(L) = 1. det(U) is the product of its diagonal
     * elements. det(P)=det(P^T) is the sgn of P. */
    for (int i = 0; i < this->n; i++)
        det *= LU[i][i];
    /* check sgn of P. write func in util. its used elsewhere also. */
    return det;
}

FMatrix FMatrix::copy() const
{
    vector<vector<GF_element>> m(this->n, vector<GF_element>(this->n));

    for (int row = 0; row < this->n; row++)
        for (int col = 0; col < this->n; col++)
            m[row][col] = this->m(row,col);

    return FMatrix(this->n, m);
}
