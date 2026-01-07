### ArcaneFEM interface to PETSc cheat  sheet

## Axl parameters

Following abbreviations can be useful for reading the table below.

- [D] - Default
- [I] - Integer type
- [R] - Real type
- [S] - String type

| Parameter | Use | Comment |
| --------- | --- | ------- |
| `solver` | Applied Krylov solver. Corresponds to the `-ksp_type` PETSc parameter. | [S] [D] = "cg"<br />**Note**: You can find all possible options on [this page](https://petsc.org/release/manual/ksp/#tab-kspdefaults). |
| `pc-type` | Applied preconditioner. Corresponds to the `-pc_type` PETSc parameter. | [S] [D] = "jacobi"<br />**Note**: You can find all possible options on [this page](https://petsc.org/release/manual/ksp/#tab-pcdefaults). |
| `rtol` | Relative convergence tolerance for the Krylov Solver. | [R] [D] = 1.0e-7 |
| `atol` | Absolute convergence tolerance for the Krylov Solver. | [R] [D]  = 0.0<br />set `rtol` to 0 if you want to use `atol` |
| `max-iter` | Maximum Krylov iterations. | [I] [D] = 1000 |

## Matrix and Vector types

The matrix and vector types use  internally by PETSc are deducted by default by a number of parameters:

- Is the program running on a single process or on multi processes ?
- Is the program running on GPU(s) (`AcceleratorRuntime` option) ?

**Notes**:
- Before using a certain type of matrix on GPU, check if it is supported on your machine. You can find this information on [this page](https://petsc.org/release/overview/gpu_roadmap).
- All possible options for matrix types can be found [here](https://petsc.org/release/overview/matrix_table).
- All possible options for vector types can be found [here](https://petsc.org/release/overview/vector_table).

## Custom options

Finally, you can specify any PETSc option you want with the `-A,petsc_flags` Arcane option. The flags specified by this option will **overwrite** the default ones (`atol`, `rtol`, `pc\_type`, `mat\_type`...). This can be useful to test some options quickly or to specify ones that are not in the `.axl` ([GAMG options](https://petsc.org/release/manualpages/PC/PCGAMG) for example).
