#include <valarray>

#include "global.hh"
#include "ematrix.hh"
#include "fmatrix.hh"
#include "gf.hh"

using namespace std;

EMatrix::EMatrix(int n, valarray<Extension_element> m): m(n, global::E.one(), m)
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
    valarray<GF_element> proj(this->n * this->n);

    for (int x = 0; x < this->n; x++)
    {
        for (int y = 0; y < this->n; y++)
            proj[x*this->n + y] = this->operator()(x,y).project();
    }

    return FMatrix(this->n, proj);
}
