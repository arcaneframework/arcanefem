## Poisson Module Benchmark Documentation

### Running the Benchmark

1.  Clone the ArcaneFEM repository to your machine.
2.  Build the Poisson module.
3.  Cd to the build folder and run the benchmark script:
```bash
$ ./run-benchmark.sh
```

### Plotting Results

The results are saved in a `results.tsv` file, which can be visualized using the `plot-results.sh` script. The script's usage is as follows:
```bash
Usage: plot-results.sh <file>
```
For example, the command
```bash
$ ./plot-results.sh benchmark-output/3D/1-mpi-instance/results.tsv
```
Will generate one plot for tests without accelerator support.

Each plot includes a title containing the parameters used for the test (based on the path of the file).

### Editing the Benchmark

The benchmark script (`arcanefem/poisson/run-benchmark.sh`) provides several global variables to customize its behavior:

#### Selecting Meshes

-   The `DIMENSIONS` variable specifies whether to run benchmarks for 2D, 3D, or both.  
    Default: `(2 3)` (both 2D and 3D).
-   Meshes are defined using `MESH_{2D|3D}` variables:
    -   **2D**: Based on `circle_cut.msh` by default.
    -   **3D**: Based on `sphere_cut.msh` by default.

You can define multiple sizes for each mesh using the `SIZES` variable. By default, 3 sizes are specified.
```bash
SIZES=("small" "medium" "large")
```

#### CPU and GPU Formats

-   `CPU_FORMATS`: An array of formats compatible with `Cpu` and `Thread` (`CT`).
-   `GPU_FORMATS`: An array of formats compatible with `Cpu`, `Thread`, and `Gpu` (`CTG`).
(The formats names should be the same as those passed by ArcaneFEM in options.)
 
The benchmark:

1.  First tests formats from `CPU_FORMATS` and `GPU_FORMATS` without the accelerator.
2.  Then reruns the `GPU_FORMATS` formats with the accelerator specified in `ACCELERATOR_RUNTIME`.

#### Configuring MPI Instances

-   Use the `MPI_N` array to define the number of MPI instances per non-accelerated test. `MPI_N_ACCELERATED` for the accelerated one.

#### Cache Warming Iterations

-   Specify the number of iterations for cache warming using the `CACHE_WARMING` variable.

### Benchmark Output

The benchmark generates results in a folder named: `benchmark-output.<CACHE_WARMING>-cw.<NUMBER>`

**Default Directory Structure:**
```bash
benchmark-output.5-cw.15
├── 2D
│   ├── 1-mpi-instance
│   │   ├── results.tsv
│   │   ├── Test.2D.1-mpi-instance.large
│   │   ├── Test.2D.1-mpi-instance.medium
│   │   └── Test.2D.1-mpi-instance.small
│   ├── 2-mpi-instance
│   │   ├── results.tsv
│   │   ├── Test.2D.2-mpi-instance.large
│   │   ├── Test.2D.2-mpi-instance.medium
│   │   └── Test.2D.2-mpi-instance.small
│   ├── cuda.1-mpi-instance
│   │   ├── results.tsv
│   │   ├── Test.2D.cuda.1-mpi-instance.large
│   │   ├── Test.2D.cuda.1-mpi-instance.medium
│   │   └── Test.2D.cuda.1-mpi-instance.small
│   └── cuda.2-mpi-instance
│       ├── results.tsv
│       ├── Test.2D.cuda.2-mpi-instance.large
│       ├── Test.2D.cuda.2-mpi-instance.medium
│       └── Test.2D.cuda.2-mpi-instance.small
├── 3D
│   ├── 1-mpi-instance
│   │   ├── results.tsv
│   │   ├── Test.3D.1-mpi-instance.large
│   │   ├── Test.3D.1-mpi-instance.medium
│   │   └── Test.3D.1-mpi-instance.small
│   ├── 2-mpi-instance
│   │   ├── results.tsv
│   │   ├── Test.3D.2-mpi-instance.large
│   │   ├── Test.3D.2-mpi-instance.medium
│   │   └── Test.3D.2-mpi-instance.small
│   ├── cuda.1-mpi-instance
│   │   ├── results.tsv
│   │   ├── Test.3D.cuda.1-mpi-instance.large
│   │   ├── Test.3D.cuda.1-mpi-instance.medium
│   │   └── Test.3D.cuda.1-mpi-instance.small
│   └── cuda.2-mpi-instance
│       ├── results.tsv
│       ├── Test.3D.cuda.2-mpi-instance.large
│       ├── Test.3D.cuda.2-mpi-instance.medium
│       └── Test.3D.cuda.2-mpi-instance.small
└── configuration.log
```
