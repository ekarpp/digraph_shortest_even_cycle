/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <valarray>
#include <vector>

#include "global.hh"
#include "fmatrix.hh"
#include "ematrix.hh"
#include "extension.hh"
#include "util.hh"

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

FMatrix FMatrix::mul_diag(GF_element e) const
{
    valarray<GF_element> m(this->n * this->n);

    for (int row = 0; row < this->n; row++)
    {
        for (int col = 0; col < this->n; col++)
        {
            if (row == col)
                m[row*this->n + col] = this->operator()(row,col) * e;
            else
                m[row*this->n + col] = this->operator()(row,col);
        }
    }

    return FMatrix(this->n, m);
}

/* PA = LU, L lower triangular with ones on diagonal
 * and U upper triangular. P permutation matrix as vector.
 * P[i] telss which row replaces row i.
 * modifies the object it is called on to L and U in single matrix. */
/* CRASHES ON SINGULAR MATRICES, maybe not? need to test */
/* modified version of the recursive leading-row-column LUP algorithm
 * presented here: https://courses.grainger.illinois.edu/cs357/fa2021/notes/ref-9-linsys.html */
vector<int> FMatrix::lup(int depth)
{
    int dim = this->n - depth;
    if (dim == 1)
        return vector<int> (this->n, depth);

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

GF_element FMatrix::det()
{
    vector<int> P = this->lup(0);
    GF_element det = global::F.one();
    /* PA = LU => det(A) = det(P^T)det(L)det(U)
     * in current implementation L has diagonal of one
     * thus det(L) = 1. det(U) is the product of its diagonal
     * elements. det(P)=det(P^T) is the sgn of P. */
    for (int i = 0; i < this->n; i++)
        det *= this->operator()(i,i);
    /* we are in characteristic two, thus -x = x and no need to check
     * sgn of P */
    return det;
}

/* uses random sampling and la grange interpolation
 * to get the polynomial determinant. */
Polynomial FMatrix::pdet(int r1, int r2) const
{
    /* determinant has deg <= 2*n - 2 */
    vector<GF_element> gamma(2*this->n - 1);
    vector<GF_element> delta(2*this->n - 1);
    uint64_t v = global::F.rem(global::randgen());

    for (int i = 0; i < 2*this->n - 1; i++)
    {
        const GF_element e(v);
        gamma[i] = e;
        /* lazy way to ensure gammas are distinct */
        v = global::F.rem(v + 2);

        FMatrix A = this->copy();
        GF_element prod = e;
        for (int col = 1; col < this->n; col++)
        {
            A.mul(r1, col, prod);
            A.mul(r2, this->n - 1 - col, prod);
            prod *= e;
        }
        delta[i] = A.det();
    }

    /* la grange */
    return util::poly_interpolation(gamma, delta);
}

FMatrix FMatrix::copy() const
{
    valarray<GF_element> m(this->n * this->n);

    for (int row = 0; row < this->n; row++)
        for (int col = 0; col < this->n; col++)
            m[row*n + col] = this->m(row,col);

    return FMatrix(this->n, m);
}
