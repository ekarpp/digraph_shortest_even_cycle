#include <vector>

#include "global.hh"
#include "fmatrix.hh"

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
