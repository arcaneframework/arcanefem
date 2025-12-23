# Introduction

**ArcaneFEM** is a collection of high-performance Finite Element Method (FEM) solvers built on the [Arcane Framework](https://github.com/arcaneframework/framework), designed for massively parallel computing on both CPU and GPU architectures. The solver showcases advanced numerical algorithms optimized for large-scale scientific computing in heterogeneous computing environments.ArcaneFEM provides a comprehensive suite of FEM-based solvers that leverage the power of modern parallel computing architectures,these solvers are engineered to handle complex engineering and scientific problems with exceptional performance and scalability. ArcaneFEM provides the following solvers in 2D and 3D:
- Laplace equation
- Poisson equation
- Fourier conduction (steady-state heat transfer)
- Electrostatics
- Potential Flow Theory (Aerodynamics)
- Transient Heat Transfer
- Bilaplacian problems
- Linear Elasticity
- Elastodynamics
- Soildynamics

### Key Features

- **Hybrid CPU-GPU Computing**: Optimized implementations for both traditional CPU clusters and modern GPU-accelerated systems
- **Massively Parallel**: Designed for large-scale distributed computing with excellent parallel efficiency
- **Domain Decomposition**: Built-in support for domain decomposition methods enabling efficient parallel execution
- **Multiple Solvers**: Includes various specialized FEM solvers for different physical phenomena (e.g., elasticity)
- **Modern Visualization**: Integrated support for ParaView visualization with HDF5 output format
- **MPI Support**: Full MPI integration for distributed memory parallelism

### Target Applications

ArcaneFEM is well-suited for:
- Heat trasfer problems
- Earthquake simulations
- Structural mechanics and elasticity problems
- Large-scale engineering simulations
- High-performance computing research

### Architecture

The solvers are built on the Arcane Framework (version 3.14.14+), which provides:
- Advanced mesh management and data structures
- Efficient memory management for large-scale problems
- Cross-platform compatibility
- Robust parallel communication layers
