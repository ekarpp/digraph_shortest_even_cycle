#include <valarray>
#include <vector>

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
 * and U upper triangular. P permutation matrix as vector.
 * P[i] telss which row replaces row i.
   modifies the object it is called on to L and U in single matrix. */
vector<int> FMatrix::lup(int depth)
{
    int dim = this->n - depth;
    if (dim == 1)
        return vector<int>(this->n, 0);

    /* bad things happen if all zero, hack and use "-1" int64_t? */
    GF_element mx = global::F.zero();
    int mxi = -1;
    for (int i = depth; i < this->n; i++)
    {
        if (this->operator()(i,depth) > mx)
        {
            mxi = i;
            mx = this->operator()(i,depth);
        }
    }

    if (mxi != depth)
        this->swap_rows(depth, mxi, depth);
    /*
       +-------+
       |M11 M12|
       |M21 M22|
       +-------+
    */

    GF_element A11_inv = this->operator()(depth, depth).inv();

    /* L21 = P22 * A21 / A11, do A21 / A11 now */
    for (int row = depth + 1; row < this->n; row++)
        this->m.mul(row, depth, A11_inv);

    /* A22 -= outerprod(A12, A22) / A11 */
    for (int row = depth + 1; row < this->n; row++)
    {
        for (int col = depth + 1; col < this->n; col++)
        {
            /* A12[col] * A21[row] */
            GF_element prod = this->operator()(depth, col)
                * this->operator()(row, depth);

            this->m.sub(row, col, prod);
        }
    }

    vector<int> P = this->lup(depth + 1);


    /* column vector of length dim - 1, starting at M[depth+1, depth] */
    valarray<GF_element> A21 = this->slice(
        depth*this->n + depth + this->n, // start
        dim - 1,                         // size
        this->n                          // stride
    );

    /* L21 do permutation */
    for (int row = depth + 1; row < this->n; row++)
        this->set(row, depth, A21[P[row] - depth - 1]);

    if (mxi == depth)
        P[depth] = depth;
    else
    {
        P[depth] = P[mxi];
        P[P[mxi]] = depth;
    }

    return P;
}

GF_element FMatrix::det() const
{
    FMatrix LU = this->copy();
    vector<int> P = LU.lup(0);
    GF_element det = global::F.one();
    /* PA = LU => det(A) = det(P^T)det(L)det(U)
     * in current implementation L has diagonal of one
     * thus det(L) = 1. det(U) is the product of its diagonal
     * elements. det(P)=det(P^T) is the sgn of P. */
    for (int i = 0; i < this->n; i++)
        det *= LU(i,i);
    /* we are in characteristic two, thus -x = x and no need to check
       sgn of P */
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
