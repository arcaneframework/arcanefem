### ArcaneFEM interface to Hypre cheat  sheet ###

| Parameter       | Type   | Use                                                  | Comment                                                      |
| --------------- | ------ | ---------------------------------------------------- | ------------------------------------------------------------ |
| `rtol`          | `Real` | relative convergence tolerance for the Krylov Solver | default = 1.0e-7                                             |
| `atol`          | `Real` | absolute convergence tolerance for the Krylov Solver | default =0.0                                                 |
| `amg-threshold` | `Real` | amg threshold strength                               | default =0.25 <br /># 2D (rep. 3D) Laplace  0.25 (rep. 0.6-0.6) is a good value<br /># elasticity, a large threshold, such as 0.9 |
| `max-iter`      | `Int`  | maximum Krylov iterations                            | default =1000                                                |
| `verbosity`     | `Int`  | verbosity for the log                                | default =2                                                   |
| `<coarsner>`    | `Int`  |                                                      | 0 : CLJP-coarsening <br />1 : classical RS coarsening (not recommended!) <br />3 : classical RS coarsening tweaked<br />6 : Falgout coarsening <br />7 : CLJP-coarsening<br />8 : PMIS-coarsening<br />10 : HMIS-coarsening<br />21 : CGC coarsening  <br />22 : CGC-E coarsening |
|                 |        |                                                      |                                                              |
|                 |        |                                                      |                                                              |
|                 |        |                                                      |                                                              |

