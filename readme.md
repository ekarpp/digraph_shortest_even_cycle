Implementation of the algorithm presented in chapter 4.2 of [The Shortest Even Cycle Problem is Tractable](https://arxiv.org/abs/2111.02992). Done as a part of Master's thesis.

# TODO
* ~n distinct random gf elements? used in fmatrix pdet and solver, current implementation is not random~ (LSFR?)
* use gaussian elimination instead of LUP for determinant
* only support GF(2^32) and GF(2^64)? (speed ups in e.g. modulus)
* AVX in GF?
* ~implement fast modulus for GF ([see here](https://dl.acm.org/doi/10.1016/j.ipl.2010.04.011))~
* further profiling after fast modulus for GF
