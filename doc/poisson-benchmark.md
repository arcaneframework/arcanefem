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
Usage: plot-results.sh <file> <mpi-n> <cache-warming> <dimension> [-cpu | -gpu]
```
For example, the command
```bash
$ ./plot-results.sh results.tsv 1 10 2D -cpu -gpu
```
will generate two plots:

-   One for tests without accelerator support.
-   One for tests with accelerator support.

Each plot includes a title containing the parameters used for the test (based on the command-line arguments).

### Editing the Benchmark

The benchmark script (`arcanefem/poisson/run-benchmark.sh`) provides several global variables to customize its behavior:

#### Selecting Meshes

-   The `DIMENSIONS` variable specifies whether to run benchmarks for 2D, 3D, or both.  
    Default: `(2 3)` (both 2D and 3D).
-   Meshes are defined using `MESH_{2D|3D}_*` variables:
    -   **2D**: Based on `circle_cut.msh` by default.
    -   **3D**: Based on `sphere_cut.msh` by default.

You can define multiple sizes for each mesh using the `SIZES` variable. By default, 3 sizes are specified.
```bash
SIZES=("small" "medium" "large")
```
Ensure that meshes are defined for all sizes and for all used dimensions.
Mesh names are dynamically resolved as `MESH` + `DIMENSION` + `SIZE` (taken from `SIZES`, in upper case), joined by `_`.

#### CPU and GPU Formats

-   `CPU_FORMATS`: An array of formats compatible with `Cpu` and `Thread` (`CT`).
-   `GPU_FORMATS`: An array of formats compatible with `Cpu`, `Thread`, and `Gpu` (`CTG`).
(The formats names should be the same as those passed by ArcaneFEM in options.)
 
The benchmark:

1.  First tests formats from `CPU_FORMATS` and `GPU_FORMATS` without the accelerator.
2.  Then reruns the `GPU_FORMATS` formats with the accelerator specified in `ACCELERATOR_RUNTIME`.

#### Configuring MPI Instances

-   Use the `CPU_CORE_NUMBERS` array to define the number of MPI instances per test.

#### Cache Warming Iterations

-   Specify the number of iterations for cache warming using the `CACHE_WARMING` variable.

### Benchmark Output

The benchmark generates results in a folder named: `benchmark-output_<CACHE_WARMING>-cw_<timestamp>`

**Default Directory Structure:**
```bash
benchmark-output_10-cw_2024-11-15_09:41:00/
├── 2D
│   ├── 1-mpi-instance
│   │   ├── results.tsv
│   │   ├── large
│   │   ├── medium
│   │   └── small
│   └── 2-mpi-instance
│       ├── results.tsv
│       ├── large
│       ├── medium
│       └── small
└── 3D
    ├── 1-mpi-instance
    │   ├── results.tsv
    │   ├── large
    │   ├── medium
    │   └── small
    └── 2-mpi-instance
        ├── results.tsv
        ├── large
        ├── medium
        └── small
```
