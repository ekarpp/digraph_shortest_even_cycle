#ifndef FMATRIX_H
#define FMATRIX_H

#include "gf.hh"
#include "matrix.hh"

class FMatrix : Matrix<GF_element>
{
    using Matrix<GF_element>::Matrix;

    void none();
};

#endif
