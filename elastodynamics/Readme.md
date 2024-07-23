# Transient elastodynamics

<img width="500" align="left" src="https://github.com/arcaneframework/arcanefem/assets/52162083/692ba9e7-5dbd-450a-ab19-e6c4a0df58a6" />

Here we deal with linear solid-mechanics governed by a system of PDE modeling the deformation of elastic bodies. The solver, here is a 2D unstructured mesh linear elasticity solver, which uses FEM to search for vector solution of displacement unknown $\mathbf{u}=(u_1,u_2)$, since transient $\mathbf{u}$ changes with time.



## Mathematics ##

#### Problem description ####

Under small elastic deformation, the steady  2D elastodynamics equation (linear elastic system) on  domain $\Omega$ reads

$$-\nabla\cdot\sigma(\mathbf{x}) +  \rho \mathbf{f}(\mathbf{x})=\rho \ddot{\mathbf{u}}(\mathbf{x}) \quad \mathbf{x}\in\Omega $$

here, $\sigma(\mathbf{x})$ is stress tensor, $\mathbf{f}(\mathbf{x})$ is the body force per unit volume, $\rho$ is the density of the material, and $\ddot{\mathbf{u}}(\mathbf{x})=\partial^2 \mathbf{u}(\mathbf{x})/\partial t^2$ is the acceleration for the displacement field $\mathbf{u}(\mathbf{x})$. The stress tensor  $\sigma(\mathbf{x})$ under isotropic elastic conditions is given by

$$ \sigma(\mathbf{x}) = \lambda(\nabla\cdot\mathbf{u}(\mathbf{x}))\mathbb{I} + \mu (\nabla\mathbf{u}(\mathbf{x}) + \left(\nabla\mathbf{u}(\mathbf{x})\right)^\text{T}) $$

here, $\lambda\in\mathbb{R}^{+}$ and $\mu\in\mathbb{R}^{+}$ are the Lame's elasticity parameters for the homogeneous material, $\mathbb{I}$ is the identity tensor, and $\mathbf{u}(\mathbf{x})$ is the displacement field vector. This governing PDE iis also knows as Navier's equation.

#### Variational formulation ####

Without entering into the details, using test function $\mathbf{v}\in\mathbb{V}$ and trial function $\mathbf{u}\in\mathbb{V}$ with $\mathbb{V}$ being the suitable FE functional space,  the varaiational formulation for the elastodynamic problem reads:

$$\int_{\Omega} \rho \ddot{\mathbf{u}}(\mathbf{x})\cdot \mathbf{v}(\mathbf{x}) +  \int_{\Omega} \lambda \nabla \cdot \mathbf{u}(\mathbf{x}) \nabla \cdot \mathbf{v}(\mathbf{x}) + 2\mu\varepsilon(\mathbf{u}(\mathbf{x})):\varepsilon(\mathbf{u}(\mathbf{x})) - \int_{\Omega}\rho\mathbf{f}(\mathbf{x})\cdot{\mathbf{v}(\mathbf{x})} - \int_{\partial\Omega_N} \mathbf{t}(\mathbf{x}) \cdot \mathbf{v}(\mathbf{x}) = 0 $$

here, $\mathbf{t}(\mathbf{x})$ is the traction vector imposed on Neumann boundary $\Omega_N$, and  $\varepsilon(\mathbf{u}(\mathbf{x})) = \varepsilon_{ij}(\mathbf{u}(\mathbf{x}))$ is the strain tensor given by

$$\varepsilon_{ij}(\mathbf{u}) = \frac{1}{2}(\frac{\partial{u}_i}{\partial{x}_j} + \frac{\partial{u}_j}{\partial{x}_i} )$$

## The code ##

This XML configuration file is used for setting up an Elastodynamics problem simulation in ArcaneFEM. Below is a detailed explanation of each section in the configuration for one such file `Test.bar.arc`.

###### Mesh Configuration

The mesh configuration section specifies the mesh file to be used in the simulation:

```xml
<meshes>
  <mesh>
    <filename>bar_dynamic.msh</filename>
  </mesh>
</meshes>
```

- **Mesh File:** Defines the mesh file (`bar_dynamic.msh`) to be used in the simulation. Note that this file should be compatible with version 4.1 `.msh` format from `Gmsh`.

###### FEM Configuration

The Finite Element Method (FEM) configuration is provided in the `Test.bar.arc`.

```xml
<fem>
  <time-discretization>Newmark-beta</time-discretization>
  <tmax>2.</tmax>
  <dt>0.08</dt>
  <rho>1.0</rho>
  <lambda>576.9230769</lambda>
  <mu>384.6153846</mu>
  <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
  <penalty>1.e64</penalty>
  <dirichlet-boundary-condition>
    <surface>surfaceleft</surface>
    <u1>0.0</u1>
    <u2>0.0</u2>
  </dirichlet-boundary-condition>
  <traction-boundary-condition>
    <surface>surfaceright</surface>
    <t2>0.01</t2>
  </traction-boundary-condition>
  <linear-system>
    <solver-backend>petsc</solver-backend>
    <preconditioner>ilu</preconditioner>
  </linear-system>
</fem>
```

Let us explain this point wise 

- **Time parameters:** The Maximum Time (tmax) time is set to `2.0`. Time Step (dt) for the simulation is set to `0.08`. Time Discretization is set via Newmark-beta

  ```xml
  <time-discretization>Newmark-beta</time-discretization>
  <tmax>2.</tmax>
  <dt>0.08</dt>
  ```

- **Material Properties:** The Density (rho) of the material is set to `1.0`.   Lame's First Parameter (lambda) is set to `576.9230769`.  The Shear Modulus (mu) is set to `384.6153846`.

  ```xml
  <rho>1.0</rho>
  <lambda>576.9230769</lambda>
  <mu>384.6153846</mu>
  ```

- **Dirichlet Boundary Condition:** Penalty  method (`Penalty`) for enforcing Dirichlet boundary conditions, with  penalty parameter for enforcing Dirichlet conditions to `1.e64`. And  the boundary condition on the specified surface (`left`) with given values for `u1` and `u2`, which we set to 0 since the end is clamped.  

  ```xml
  <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
  <penalty>1.e64</penalty>
  <dirichlet-boundary-condition>
    <surface>left</surface>
    <u1>0.0</u1>
    <u2>0.0</u2>
  </dirichlet-boundary-condition>
  ```

- **Traction Boundary Condition:** Defines the traction boundary condition on the specified surface (`surfaceright`) with a given value for `t2`.

  ```xml
  <traction-boundary-condition>
    <surface>surfaceright</surface>
    <t2>0.01</t2>
  </traction-boundary-condition>
  ```

- **Linear System Configuration:** Specifies the linear system settings, including the solver backend (`petsc`) and the preconditioner (`ilu`).

  ```xml
  <linear-system>
    <solver-backend>petsc</solver-backend>
    <preconditioner>ilu</preconditioner>
  </linear-system>
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



#### Post Process ####

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).



## Tests present in the module ##

| Test           | Arc file                          | Comment                                                      |
| -------------- | --------------------------------- | ------------------------------------------------------------ |
| elastodynamics | `Test.Elastodynamics.arc`         | Clamped bar being pulled at other end via transient load.<br />**Mesh** `dar_dynamic.msh` . Time-dicretization - **Newmark**-$\beta$<br />**No damping** is present. **Penalty** method for applying Dirichlet |
| elastodynamics | `Test.Elastodynamics.Galpha.arc`  | Clamped bar being pulled at other end via transient load.<br />**Mesh** `dar_dynamic.msh` . Time-dicretization - **Generalized**-$\alpha$<br />**No damping** is present. **RowColumnElimination** method for applying Dirichlet |
| elastodynamics | `Test.Elastodynamics.damping.arc` | Clamped bar being pulled at other end via transient load.<br />**Mesh** `dar_dynamic.msh` . Time-dicretization - **Newmark**-$\beta$<br />**Damping** is present. **Penalty** method for applying Dirichlet |
| elastodynamics | `Test.Elastodynamics.pointBC.arc` | Semi-circular section of soil, loaded via **point source** on top.<br />**Mesh** `semi-circle.msh` . Time-dicretization - **Newmark**-$\beta$<br />**No damping** is present. **RowColumnElimination** method for applying Dirichlet |
|                |                                   |                                                              |

