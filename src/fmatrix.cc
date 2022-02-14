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
