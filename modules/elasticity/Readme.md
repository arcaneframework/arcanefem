# linear elasticity

<img width="500" align="left" src="https://github.com/arcaneframework/arcanefem/assets/52162083/eb970ece-5fd3-4862-9b93-e8930a103ae9" />

Here, we focus on linear solid mechanics, which involves analyzing the behavior of elastic bodies through a system of partial differential equations (PDEs) that govern their deformation. The solver at hand is specifically designed for 3D/2D unstructured mesh problems and operates as a linear elasticity solver utilizing the finite element method (FEM).

The primary objective of the solver is to determine the vector solution for the unknown displacements within the elastic body. By utilizing FEM, it efficiently searches for a solution that satisfies the given boundary conditions and accurately captures the deformation characteristics of the system. This solver serves as a valuable tool for investigating and understanding the mechanical behavior of elastic materials in various engineering applications.

## Mathematics

#### Problem description

Under small elastic deformation, the steady elastic deformation equation (linear elastic system) on  domain $\Omega$ reads

$$
-\nabla\cdot\sigma(\mathbf{x})=\mathbf{f}(\mathbf{x}) \quad \mathbf{x}\in\Omega
$$

here, $\sigma(\mathbf{x})$ is stress tensor, $\mathbf{f}(\mathbf{x})$ is the body force per unit volume. The stress tensor  $\sigma(\mathbf{x})$ under isotropic elastic conditions is given by

$$
\sigma(\mathbf{x}) = \lambda(\nabla\cdot\mathbf{u}(\mathbf{x}))\mathbb{I} + \mu (\nabla\mathbf{u}(\mathbf{x}) + \left(\nabla\mathbf{u}(\mathbf{x})\right)^\text{T})
$$

here, $\lambda\in\mathbb{R}^{+}$ and $\mu\in\mathbb{R}^{+}$ are the Lame's elasticity parameters for the homogeneous material, $\mathbb{I}$ is the identity tensor, and $\mathbf{u}(\mathbf{x})$ is the displacement field vector. This governing PDE iis also knows as Navier's equation.

#### Variational formulation

Without entering into the details, the variational formulation for the Navier's equation reads

$$
\int_{\Omega} \lambda \nabla \cdot \mathbf{u}(\mathbf{x}) \nabla \cdot \mathbf{v}(\mathbf{x}) + 2\mu\varepsilon(\mathbf{u}(\mathbf{x})):\varepsilon(\mathbf{u}(\mathbf{x})) - \int_{\Omega}\mathbf{f}(\mathbf{x})\cdot{\mathbf{v}(\mathbf{x})} - \int_{\partial\Omega_N} \mathbf{t}(\mathbf{x}) \cdot \mathbf{v}(\mathbf{x}) = 0
$$

here, $\mathbf{t}(\mathbf{x})$ is the traction vector imposed on Neumann boundary $\Omega_N$, and  $\varepsilon(\mathbf{u}(\mathbf{x})) = \varepsilon_{ij}(\mathbf{u}(\mathbf{x}))$ is the strain tensor given by

$$
\varepsilon_{ij}(\mathbf{u}) = \frac{1}{2}(\frac{\partial{u}_i}{\partial{x}_j} + \frac{\partial{u}_j}{\partial{x}_i} )
$$

## Running the solver

Here we will provide a simple 2D bar problem example. This XML configuration file is used for setting up an Elastodynamics problem simulation in ArcaneFEM. Below is a detailed explanation of each section in the configuration for one of the tests `Test.bar.arc`.

###### Mesh Configuration

The mesh configuration section specifies the mesh file to be used in the simulation:

```xml
<meshes>
  <mesh>
    <filename>meshes/bar.msh</filename>
      <subdivider>
        <nb-subdivision>0</nb-subdivision>
      </subdivider>
  </mesh>
</meshes>

```

- **Mesh File:** Defines the mesh file (`bar.msh`) to be used in the simulation. Note that this file should be compatible with version 4.1 `.msh` format from `Gmsh`.
- **Subdivision**: `nb-subdivision` here allows to split the mesh to produce a finer mesh from a coarse mesh. It is practical to use for producing meshes for large-scale simulations.

###### FEM Configuration

The Finite Element Method (FEM) configuration is provided in the `Test.bar.arc`.

```xml
<fem>
  <result-file>test_elasticity_results.txt</result-file>
  <E>21.0e5</E>
  <nu>0.28</nu>
  <f>NULL -1.0</f>
  <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
  <dirichlet-boundary-condition>
    <surface>left</surface>
    <u>0.0 0.0</u>
  </dirichlet-boundary-condition>
</fem>
```

Let us explain this point wise

- **Result File:** Specifies the file for validation test (OPTIONAL). Use only if you want to compare solutions against known solutions

  ```xml
  <result-file>test_elasticity_results.txt</result-file>
  ```
- **Material Properties:** The Young's Modulus (E) for the material, defined as `21.0e5`.  The Poisson's Ratio (nu) for the material, defined as `0.28`.

  ```xml
  <E>21.0e5</E>
  <nu>0.28</nu>
  ```
- **Source Term:** The source term or dody force is present in in Y direction hence (f2) in the PDE, set to `-1.0`.

  ```xml
  <f>NULL -1.0</f>
  ```
- **Dirichlet Method:** Specifies the method (`Penalty`) for enforcing Dirichlet boundary conditions. And  the boundary condition on the specified surface (`left`) with given values for `u1` and `u2`, which we set to 0 since the end is clamped.

  ```xml
  <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
  <dirichlet-boundary-condition>
    <surface>left</surface>
    <u>0.0 0.0</u>
  </dirichlet-boundary-condition>
  ```

###### Post-Processing Configuration

The post-processing configuration is specified to control how and when results are saved:

```xml
<arcane-post-processing>
  <output-period>1</output-period>
  <format name="VtkHdfV2PostProcessor" />
  <output>
    <variable>U</variable>
  </output>
</arcane-post-processing>
```

- **Output Period:** Sets the interval at which results are saved.
- **Format:** Specifies the format for the output files (`VtkHdfV2PostProcessor`).
- **Output Variables:** Lists the variables (`U`) which is the displacement vector to be included in the output.

#### Post Process

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).

#### Tests available in this module

The tests are present in the form of `.arc` files with a prefix `Test.`:

| Name                                 | Dimension | Boundary Condition                                          | Solver          | Comment                                                       |
| ------------------------------------ | --------- | ----------------------------------------------------------- | --------------- | ------------------------------------------------------------- |
| bar                                  | 2D        | Clamped Dirichlet + Null flux<br />Body force (Source Term) | Default (PETSc) | - Serves as validation test                                   |
| bar.DirichletViaRowElimination       | 2D        | Clamped Dirichlet + Null flux<br />Body force (Source Term) | PETSc           | - Row Elimination for BC<br />- GMRES with ILU(0) Solver      |
| bar.DirichletViaRowColumnElimination | 2D        | Clamped Dirichlet + Null flux<br />Body force (Source Term) | PETSc           | - Row and Column Elimination for BC<br />- CG with AMG Solver |
| bar.traction                         | 2D        | Clamped Dirichlet<br />Traction (Neumann)                   | Default (PETSc) |                                                               |
|                                      |           |                                                             |                 |                                                               |
