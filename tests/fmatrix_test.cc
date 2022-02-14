#include <iostream>
#include <vector>

#include "fmatrix_test.hh"
#include "../src/fmatrix.hh"
#include "../src/global.hh"
#include "../src/gf.hh"

using namespace std;

FMatrix_test::FMatrix_test(int dim)
{
    this->dim = dim;

    cout << "---------------" << endl;
    cout << "TESTING FMATRIX" << endl;
    cout << "---------------" << endl << endl;

    this->run();
}

FMatrix FMatrix_test::random_matrix()
{
    vector<vector<GF_element>> m(this->dim, vector<GF_element>(this->dim));
    for (int i = 0; i < this->dim; i++)
        for (int j = 0; j < this->dim; j++)
            m[i][j] = global::F.random();
    return FMatrix(this->dim, m);
}

void FMatrix_test::test_addition()
{
    cout << "addition" << endl;
    int err = 0;
    for (int i = 0; i < this->tests; i++)
    {
        FMatrix a = this->random_matrix();
        FMatrix b = this->random_matrix();

        GF_element sum = global::F.zero();
        for (int i = 0; i < this->dim; i++)
            for (int j = 0; j < this->dim; j++)
                sum = sum + a[i][j] + b[i][j];
        FMatrix c = a + b;
        //if (sum != a + b)
        //    err++;
    }
    end_test(err);
}

void FMatrix_test::test_multiplication()
{

}
