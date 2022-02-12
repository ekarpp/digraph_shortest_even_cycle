#include "../src/global.hh"
#include "gf_test.hh"

util::rand64bit global::randgen;

int main(void)
{
    global::randgen.init(10);
    GF_test t(10);

    t.test_add_inverse();
    t.test_mul_commutative();
    t.test_mul_inverse();
    t.test_mul_id();
    return 0;
}
