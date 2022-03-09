The code of [Low-Rank Solution of Lyapunov Equations](https://doi.org/10.1137/S0895479801384937)

**Algorithm 1.** Alternating direction implicit algorithm.

ADI.m

**Algorithm 2.** The CF–ADI algorithm.

CFADI.m



The translation of the paper is in [Low-Rank Solution of Lyapunov Equations（一）ADI算法](https://blog.csdn.net/qq_34179307/article/details/123149380?spm=1001.2014.3001.5501)、[（二）CF-ADI算法](https://blog.csdn.net/qq_34179307/article/details/123227974)



Input the values of A and B, referring to the content of [A cyclic low-rank Smith method for large sparse Lyapunov equations](https://doi.org/10.1137/S1064827598347666), using the simplest form of A matrix -- symmetric triangular matrix.



We didn't write the code for m<1, it was too complicated.