# Solving Laplace equation with ArcaneFEM #
<img align="left" width="400" src="https://github.com/arcaneframework/arcanefem/assets/52162083/be3d2ea6-bfb7-42d9-b82e-a62509a498f8" />


Here, we employ Arcane to solve the Laplace equation, which stands as one of the fundamental partial differential equations (PDEs). The provided code presents a straightforward implementation of a 2D unstructured mesh Galerkin finite element method (FEM) solver.

The Laplace equation, also known as the harmonic equation, appears in various scientific and engineering applications. It plays a crucial role in fields such as physics, engineering, and mathematics. By solving the Laplace equation, we can gain insights into phenomena such as steady-state heat conduction, electrostatics, fluid flow, potential fields, and more. The code showcased here provides a foundation for tackling problems governed by this fundamental equation.


#### Problem description ####

The 2D Laplace equation is solved for a closed meshed domain $\Omega^h$ in order to know the Laplace solution $u(x,y)$ within the domain. The equation reads

$$\frac{\partial}{\partial x}\left(\frac{\partial u}{\partial x} \right) + \frac{\partial}{\partial y}\left(\frac{\partial u}{\partial y} \right) = 0  \quad \forall (x,y)\in\Omega^h $$

or in a more compact forms

$$\nabla(\nabla u)= 0 \quad \forall (x,y)\in\Omega^h.$$ or

$$\nabla^2 u= 0 \quad \forall (x,y)\in\Omega^h.$$ or

$$\Delta^2 u= 0 \quad \forall (x,y)\in\Omega^h.$$ or

To complete the problem description,  three first type (Dirichlet) boundary conditions are applied to this problem:

$u = 50.0 \quad \forall(x,y)\in\partial\Omega^h_{\text{inner}}\subset\partial \Omega^h,$

$u = 20.0 \quad \forall(x,y)\in\partial\Omega^h_{\text{outer}}\subset\partial \Omega^h,$ and

We work with approximation, $\lambda$ is homogeneous $\lambda : \Omega^h \in \mathbb{R}^{+}$, in this case  the variational formulation in $H^1_{0}(\Omega) \subset H^1{\Omega}$  reads

search FEM trial function $u^h(x,y)$ satisfying

$$- \int_{\Omega^h}\lambda\nabla u^h \nabla  v^h + \int_{\partial\Omega_N} (\overline{q} \cdot \mathbf{n}) v^h v^h = 0 \quad \forall v^h\in H^1_0(\Omega^h)$$

given

$u^h=50.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{inner}}$,

$u^h=20.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{outer}}$,

$\int_{\Omega^h_{{N}}}(\mathbf{q} \cdot \mathbf{n}) v^h=0$

## The code ##

#### Mesh ####

The mesh `plancher.msh` is provided in the `Test.ring.arc` file

```xml
  <meshes>
    <mesh>
      <filename>ring.msh</filename>
    </mesh>
  </meshes>
```

Please not that use version 4.1 `.msh` file from `Gmsh`.

#### Boundary conditions ####

The Dirichlet (constant $u$) boundary conditions  are provided in `Test.ring.arc` file

```xml
    <dirichlet-boundary-condition>
      <surface>inner</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>outer</surface>
      <value>20.0</value>
    </dirichlet-boundary-condition>
```

So in the snippet above, three Dirichlet conditions are applied ($50, 20.0$)  on three borders ('inner', 'outer') for the loaded mes `ring.msh`.

The Neumann  boundary conditions  are absent but could be provided in such a way

```xml
    <neumann-boundary-condition>
      <surface>outer</surface>
      <value>16.0</value>
    </neumann-boundary-condition>
```



#### Post Process ####

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).



#### Tests available in this module ####

The tests are present in the form of `.arc` files with a prefix `Test.`:

| Name                   | Dimension | Boundary Condition          | Solver               | Comment                     |
| ---------------------- | --------- | --------------------------- | -------------------- | --------------------------- |
| L-shape                | 3D        | Dirichlet + Null flux       | Default (PETSc)      | - Serves as validation test |
| ring                   | 2D        | Dirichlet only              | Default (PETSc)      |                             |
| PointDirichlet         | 2D        | Point Dirichlet + Null flux | Sequential Direct LU | - Serves as validation test |
| PointDirichlet-refined | 2D        | Point Dirichlet + Null flux | Default (PETSc)      | - Refined version of test   |
|                        |           |                             |                      |                             |

####  Point loading example ####

<img src="https://github.com/arcaneframework/arcanefem/assets/52162083/979b2bd6-9e54-4fae-b6d7-f2b5a232b4dc" alt="Test_2_large_new" style="zoom:67%;" />

