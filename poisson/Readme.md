# Solving Poisson equation with ArcaneFEM #

<img align="left" width="400" src="https://github.com/arcaneframework/arcanefem/assets/52162083/3646fbb9-5c82-4807-b025-7f3b4c899ca7" alt="poisson_1_large" />

Here, we utilize ArcaneFEM to solve the Poisson equation, which is a fundamental elliptic partial differential equation (PDE). The provided code demonstrates a straightforward implementation of a 2D/3D unstructured mesh Galerkin finite element method (FEM) solver on an L-shaped domain. Although we shall explain here only 2D for keeping the text simple.

The Poisson equation is encountered in various physical scenarios, including heat conduction, substance diffusion, membrane elasticity, inviscid fluid flow, electrostatics, twisting of elastic rods, and water waves. It serves as a vital tool for modeling and understanding these phenomena.

## Theory ##

#### Problem description ####

The 2D Poisson's equation is solved for a closed meshed domain $\Omega^h$ in order to know the solution $u(x,y)$ within the domain. The equation reads

$$\frac{\partial}{\partial x}\left( \frac{\partial u}{\partial x} \right) + \frac{\partial}{\partial y}\left( \frac{\partial u}{\partial y} \right) = {\mathcal{f}}   \quad \forall (x,y)\in\Omega^h $$

or in a more compact form

$$\nabla^2 u = {\mathcal{f}} \quad \forall (x,y)\in\Omega^h.$$



To complete the problem description,   first type (Dirichlet) boundary conditions is applied to this problem:

$u = 0.0 \quad \forall(x,y)\in\partial\Omega^h_{\text{boundary}}\subset\partial \Omega^h,$

Finally the right hand side source is present within the domain

${\mathcal{f}}=-1$



In this case  the FEM variational formulation in $H^1_{0}(\Omega) \subset H^1{\Omega}$  reads

search FEM trial function $u^h(x,y)$ satisfying

$$- \int_{\Omega^h}\nabla u^h \nabla  v^h + \int_{\partial\Omega_N} (\overline{q} \cdot \mathbf{n}) v^h + \int_{\Omega^h}{\mathcal{f}} v^h = 0 \quad \forall v^h\in H^1_0(\Omega^h)$$

given

$u^h=0.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{boundary}}$,

$\int_{\Omega^h_{\text{N}}}(\mathbf{q} \cdot \mathbf{n}) v^h=0$ since no Neumann BC is present,

$\int_{\Omega^h}{\mathcal{f}} v^h=1\times10^5$, and

## The code ##

#### properties ###

The value of constant source term $\mathcal{f}$  can be provided in  `Test.L-shape.2D.arc` file

```xml
  <fem>
    <f>-1</f>
  </fem>
```

#### Mesh ####

The mesh `L-shape.msh` is provided in the `Test.L-shape.2D.arc` file

```xml
  <meshes>
    <mesh>
      <filename>L-shape.msh</filename>
    </mesh>
  </meshes>
```

Note, here `L-shape.msh` is a 2D mesh, if any other 3D mesh was loaded ArcaneFEM will run 3D calculations it is as simple as that.

Please not that use version 4.1 `.msh` file from `Gmsh`.

#### Boundary conditions ####

The Dirichlet boundary conditions  are provided in `Test.L-shape.2D.arc` file

```xml
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
```

So in the snippet above, three Dirichlet condition $u=0$ is  applied to border ('boundary') which is a group of edges in the mesh file `L-shape.msh`.

If needed, the Neumann  boundary conditions  can also be provided in `Test.L-shape.2D.arc` file

```xml
    <neumann-boundary-condition>
      <surface>Neumann</surface>
      <value>0.0</value>
    </neumann-boundary-condition>
```



#### Post Process ####

<img align="left" width="400" src="https://github.com/arcaneframework/arcanefem/assets/52162083/a8d114e1-5589-4efd-88fd-84b398acab84" alt="poisson_1_large" />

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).





#### Tests available in this module ####

The tests are present in the form of `.arc` files with a prefix `Test.`:

| Name          | Dimension | Boundary Condition                                   | Solver               | Comment                                                     |
| ------------- | --------- | ---------------------------------------------------- | -------------------- | ----------------------------------------------------------- |
| L-shape.2D    | 2D        | Homogeneous Dirichlet <br />Homogeneous source term  | Default (PETSc)      | - Serves as validation test                                 |
| L-shape.3D    | 3D        | Dirichlet + Null flux <br />Homogeneous source term  | PETSc                | - Serves as validation test<br />- Uses BLCSR matrix format |
| sphere.3D     | 3D        | Dirichlet only<br />Homogeneous source term          | HYPRE (ArcaneFEM)    | - Uses BLCSR matrix format                                  |
| direct-solver | 2D        | Homogeneous Dirichlet <br />Homogeneous source term  | Sequential Direct LU |                                                             |
| neumann       | 2D        | Neumann only                                         | Default (PETSc)      | - Serves as validation test                                 |
| porous        | 2D        | Multiple Dirichlet only<br />Homogeneous source term | PETSc                | - Used for Benchmarking                                     |
| trilinos      | 2D        | Homogeneous Dirichlet <br />Homogeneous source term  | TRILINOS             | - Used to test TRILINOS                                     |
| petsc         | 2D        | Homogeneous Dirichlet <br />Homogeneous source term  | PETSc                | - Serves as validation test                                 |
| hypre         | 2D        | Homogeneous Dirichlet <br />Homogeneous source term  | HYPRE (Arcane)       | - Serves as validation test                                 |
| hypre_direct  | 2D        | Homogeneous Dirichlet <br />Homogeneous source term  | HYPRE (ArcaneFEM)    | - Serves as validation test<br />- Uses BLCSR matrix format |
|               |           |                                                      |                      |                                                             |



### Time analysis ###
By setting the REGISTER_TIME flag to ON during the compilation, it is possible to generate a `timer.txt` file during the execution which contains the execution time of the
different parts of poisson.

Here is an example of compilation with this flag :
~~~{sh}
ARCANE_INSTALL_DIR=/path/to/arcane/installation
BUILD_DIR=/tmp/build
SOURCE_PATH=/path/to/sources
cmake -S ${SOURCE_PATH} -B ${BUILD_DIR} -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_DIR} -DREGISTER_TIME=ON
cmake --build ${BUILD_DIR}
~~~

