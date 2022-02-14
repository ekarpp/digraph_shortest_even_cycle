#include <vector>

#include "global.hh"
#include "ematrix.hh"
#include "fmatrix.hh"
#include "gf.hh"

using namespace std;

EMatrix::EMatrix(int n, vector<vector<Extension_element>> m): m(n, global::E.one(), m)
{
    this->n = n;
}

EMatrix::EMatrix(Matrix<Extension_element> m)
{
    this->n = m.get_n();
    this->m = m;
}

EMatrix EMatrix::operator+(const EMatrix &other) const
{
    return EMatrix(this->m + other.get_m());
}

EMatrix EMatrix::operator-(const EMatrix &other) const
{
    return EMatrix(this->m - other.get_m());
}

EMatrix EMatrix::operator*(const EMatrix &other) const
{
    return EMatrix(this->m * other.get_m());
}

FMatrix EMatrix::project() const
{
    vector<vector<GF_element>> proj(this->n, vector<GF_element>(this->n));

    for (int x = 0; x < this->n; x++)
    {
        const vector<Extension_element> row = this->m[x];
        for (int y = 0; y < this->n; y++)
            proj[x][y] = row[y].project();
    }

    return FMatrix(this->n, proj);
}
