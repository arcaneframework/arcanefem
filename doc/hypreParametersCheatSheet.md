### ArcaneFEM interface to Hypre cheat  sheet ###

Following abbreviations can be useful for reading the table below.

- [D] - Default
- [I] - Integer type
- [R] - Real type

| Parameter         | Use                                                  | Comment                                                      |
| ----------------- | ---------------------------------------------------- | ------------------------------------------------------------ |
| `rtol`            | relative convergence tolerance for the Krylov Solver | [R] [D] = 1.0e-7                                             |
| `atol`            | absolute convergence tolerance for the Krylov Solver | [R] [D]  = 0.0<br />set `rtol` to 0 if you want to use `arol` |
| `max-iter`        | maximum Krylov iterations                            | [I] [D] = 1000                                               |
| `verbosity`       | verbosity for the log                                | [I] [D] = 2                                                  |
| `amg-coarsener`   | amg parallel coarsening algorithm                    | [I] [D] = 8<br />0   : CLJP <br />3   : classical RS <br />6   : Falgout  <br />8   : PMIS<br />10 : HMIS<br />21 : CGC   <br />22 : CGC-E<br />**Note** GPU supported choice 8 |
| `amg-threshold`   | amg threshold strength                               | [R] [D] =0.25 <br /># 2D (rep. 3D) Laplace  0.25 (rep. 0.6-0.6) is a good value<br /># Elasticity, a large threshold, such as 0.9 |
| `amg-interp-type` | amg interpolation operator                           | [I] [D] = 6<br />0   : classical <br />1   : LS  (for use with GSMG) <br />2   : classical (hyperbolic PDEs)<br />3   : direct with separation of weights <br />4   : multipass<br />5   : multipass with separation of weights<br />6   : extended+i<br />7   : extended+i no common C neighbor<br />8   : standard <br />9   : standard with separation of weights<br />10 : classical block (use for nodal systems)<br />11  : 10 with diagonalized diagonal blocks  <br />12  : FF interpolation<br />13  : FF1 interpolation<br />14  : extended<br />15  : adaptive weights<br />16  : extended  in matrix-matrix form<br />17  : extended+i in matrix-matrix form<br />18  : extended+e in matrix-matrix form<br />**Note** GPU supported choice 3, 6, 14, 15, 18 |
| `amg-smoother`    | amg smoother to be used                              | [I] [D] = 6<br />0 : Jacobi<br/>3 : hybrid GS or SOR, forward solve<br/>4 : hybrid GS or SOR, backward solve<br/>6 : hybrid symmetric GS or SSOR<br/>7 : Jacobi with Matvec<br />8 : $l_1$-scaled hybrid symmetric GS<br/>9 : Gaussian elim. (on coarsest level)<br/>10 : On-proc. direct forward solve for matrices with triangular structure<br/>11 : Two Stage approx. to GS. Uses lower part diagonal matrix<br/>12 : 11 and a second iteration for error approx.<br/>13 : $l_1$ Gauss-Seidel, forward solve<br/>14 : $l_1$ Gauss-Seidel, backward solve<br/>16 : Chebyshev<br/>17 : FCF-Jacobi<br/>18 : $l_1$-scaled jacobi<br/>30 : Kaczmarz<br/>88: 8 with a convergent l1-term<br/>89: Symmetric l1-hybrid GS (i.e., 13 followed by 14)<br/>98 : LU with pivoting<br />**Note** GPU supported choice 3, 4, 6, 7, 18, 11, 12 |
|                   |                                                      |                                                              |

