# Solving Poisson equation with FEM and Arcane #

Here Poisson equation, which is one of the basics elliptic PDEs, is solved using FEM in Arcane. The code here is a simple 2D unstructured mesh Galerkin FEM solver. The Poisson equation arises in numerous physical contexts, e.g., heat conduction, diffusion of substances, membrane elasticity,  inviscid fluid flow, electrostatics, twisting of elastic rods, and water waves.


![poisson_1_small](https://github.com/arcaneframework/arcanefem/assets/52162083/2864be45-1e21-47be-90df-db1fbfc7bee8)

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

The value of constant source term $\mathcal{f}$  can be provided in  `Test.poisson.arc` file

```xml
  <fem>
    <f>-1</f>
  </fem>
```

#### Mesh ####

The mesh `L-shape.msh` is provided in the `Test.poisson.arc` file

```xml
  <meshes>
    <mesh>
      <filename>L-shape.msh</filename>
    </mesh>
  </meshes>
```

Please not that use version 4.1 `.msh` file from `Gmsh`.

#### Boundary conditions ####

The Dirichlet boundary conditions  are provided in `Test.poisson.arc` file

```xml
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
```

So in the snippet above, three Dirichlet condition $u=0$ is  applied to border ('boundary') which is a group of edges in the mesh file `L-shape.msh`.

If needed, the Neumann  boundary conditions  can also be provided in `Test.poission.arc` file

```xml
    <neumann-boundary-condition>
      <surface>Neumann</surface>
      <value>0.0</value>
    </neumann-boundary-condition>
```



#### Post Process ####

For post processing the `ensight.case` file is outputted, which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).
