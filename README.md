# ArcaneFEM

**High-Performance Finite Element Method Solvers with CPU/GPU Parallelism**

ArcaneFEM provides production-ready Finite Element Method (FEM) solvers built on the [Arcane Framework](https://github.com/arcaneframework/framework). Designed for modern HPC environments, these solvers deliver optimized performance across diverse parallel computing architectures: multi-CPU, multi-GPU, and hybrid CPU-GPU configurations.

##### Documentation & Resources

- **[ArcaneFEM Documentation](https://arcaneframework.github.io/arcanefem/)** - Still in progress but usable user guide
- **[GitHub Repository](https://github.com/arcaneframework/arcanefem)** - Source code and issue tracking
- **[Arcane Framework Docs](https://arcaneframework.github.io/arcane/userdoc/html/index.html)** - Core framework reference

##### Key Features

- **Flexible Parallelism**: Seamlessly run on CPUs, GPUs, or heterogeneous CPU-GPU systems
- **Multiple Physics Modules**: Includes solvers for elasticity, heat transfer, and more
- **Production-Ready**: Optimized algorithms for real-world engineering simulations
- **Modern Visualization**: Native support for ParaView via VTKHDF5 and Ensight formats

------

## Installation Notes

[Detailed Installation procedure](https://arcaneframework.github.io/arcanefem/#/install)

#### Required Dependencies

- **[Arcane Framework](https://arcaneframework.github.io/arcane/userdoc/html/d7/d94/arcanedoc_build_install.html)** - Core computational framework
- **Linear Solver Library** (at least one):
  - **HYPRE** (recommended for CPU and GPU parallelism)
  - **PETSc** (recommended for CPU and GPU parallelism)
  - **Trilinos**

Refer to the [Arcane Installation Guide](https://arcaneframework.github.io/arcane/userdoc/html/d7/d94/arcanedoc_build_install.html) for detailed compilation instructions. **Important**: Configure Arcane with HYPRE, PETSc, or Trilinos support to unlock ArcaneFEM's full capabilities.

#### Building ArcaneFEM

Assuming Arcane Framework is already installed:

```bash
# Configure paths
export ARCANE_INSTALL_DIR=/path/to/arcane/installation
export ARCANEFEM_INSTALL_DIR=${HOME}/ArcaneFEM/install
export SOURCE_PATH=/path/to/ArcaneFEM/sources

# Configure build with CMake
cmake -S ${SOURCE_PATH} \
      -B ${ARCANEFEM_INSTALL_DIR} \
      -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_DIR}

# Compile and install
cmake --build ${ARCANEFEM_INSTALL_DIR}
```

------

## Quick Start Guide

#### Running Your First Simulation

After compilation, navigate to a solver module. For example, the elasticity solver:

```bash
cd ${ARCANEFEM_INSTALL_DIR}/modules/elasticity
```

> **Tip**: Explore other physics modules in `${ARCANEFEM_INSTALL_DIR}/modules/`

Each module includes example input files in its `inputs/` directory.

#### Execution Modes

**Sequential (Single Core)**

```bash
./Elasticity ./inputs/Test.Elasticity.arc
```

**Parallel CPU (Domain Decomposition)**

```bash
# Run on 4 CPU cores
mpirun -n 4 ./Elasticity ./inputs/Test.Elasticity.arc
```

**GPU Accelerated**

```bash
# Single NVIDIA GPU
mpirun -n 1 -A,AcceleratorRuntime=cuda ./Elasticity ./inputs/Test.Elasticity.arc

# Single AMD GPU
mpirun -n 1 -A,AcceleratorRuntime=hip ./Elasticity ./inputs/Test.Elasticity.arc
```

> **Note**: Replace `1` with the number of available GPUs

**Hybrid CPU-GPU**

```bash
# 8 CPU cores + 1 GPU
mpirun -n 8 -A,AcceleratorRuntime=cuda ./Elasticity ./inputs/Test.Elasticity.arc
```

For advanced runtime options, consult the [Arcane Launcher Documentation](https://arcaneframework.github.io/arcane/userdoc/html/d8/dd6/arcanedoc_execution_launcher.html).

#### Visualization

ArcaneFEM outputs results in modern visualization formats:

- **VTKHDF5** (`*.hdf`) - Recommended for large datasets
- **Ensight** (`*.case`) - Legacy format support

Results are written to the `output/` directory within each module.

##### Viewing Results with ParaView

```bash
# Open VTKHDF5 output
paraview ${ARCANEFEM_INSTALL_DIR}/modules/elasticity/output/depouillement/vtkhdfv2/Mesh0.hdf
```

**Requirements**:

- ParaView ≥ 5.12
- Arcane compiled with MPI support for HDF5

------

## Gallery

*Coming soon: Simulation examples, benchmark results, and application showcases*

------

## Getting Help

- **Issues**: Report bugs or request features on [GitHub Issues](https://github.com/arcaneframework/arcanefem/issues)
- **Documentation**: Detailed solver guides at [arcaneframework.github.io/arcanefem](https://arcaneframework.github.io/arcanefem/)
