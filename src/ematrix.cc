#include <vector>

#include "global.hh"
#include "ematrix.hh"

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
