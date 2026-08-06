# Solving Poisson equation with Discontinuous Galerkin method ArcaneFEM

<img align="left" width="400" src="https://github.com/arcaneframework/arcanefem/assets/52162083/3646fbb9-5c82-4807-b025-7f3b4c899ca7" alt="poisson_1_large" />

Here, we utilize ArcaneFEM to solve the Poisson equation, which is a fundamental elliptic partial differential equation (PDE). The provided code demonstrates a straightforward implementation of a 2D/3D unstructured mesh Discontinous Galerkin finite element method (DG-FEM) solver on arbitary meshes. Although we shall explain here only 2D for keeping the text simple. 

The Poisson equation is encountered in various physical scenarios, including heat conduction, substance diffusion, membrane elasticity, inviscid fluid flow, electrostatics, twisting of elastic rods, and water waves. It serves as a vital tool for modeling and understanding these phenomena.

## Theory

#### Problem description

The 2D Poisson's equation is solved for a closed meshed domain $\Omega^h$ in order to know the solution $u(x,y)$ within the domain. The equation reads

$$
\frac{\partial}{\partial x}\left( \frac{\partial u}{\partial x} \right) + \frac{\partial}{\partial y}\left( \frac{\partial u}{\partial y} \right) = {\mathcal{f}}   \quad \forall (x,y)\in\Omega^h
$$

or in a more compact form

$$
\nabla^2 u = {\mathcal{f}} \quad \forall (x,y)\in\Omega^h.
$$

To complete the problem description,   first type (Dirichlet) boundary conditions is applied to this problem:

$u = 0.5 \quad \forall(x,y)\in\partial\Omega^h_{\text{horizontal}}\subset\partial \Omega^h,$

Finally the right hand side source is present within the domain

${\mathcal{f}}=5.5$

In this case the DG-FEM variational formulation in $H^1(\Omega)$ reads

search DG-FEM trial function $u^h(x,y)$ satisfying

$$ \sum_K \int_{K} \nabla u^h \nabla v^h - \sum_F \int_{F} (\{\nabla u^h \cdot \mathbf{n}\} [v^h] + \{\nabla v^h \cdot \mathbf{n}\} [u^h] - \alpha h_F^{-1} [u^h][v^h]) + \sum_{F\subset\partial\Omega} \int_F g v^h = \int_{\Omega} {\mathcal{f}} v^h \quad \forall v^h \in H^1(\Omega)$$

given

$u^h = g \quad \forall (x,y)\in\partial\Omega^h_{\text{Dirichlet}}\subset\partial \Omega^h,$

where $K$ are elements, $F$ faces, $\{\cdot\}$ average, $[\cdot]$ jump, $\alpha$ penalty parameter, $h_F$ face size.

## The code ##

Once you compile the Poisson_dg module. Look out for three things in the folder: 
- `Poisson_dg` this is the executable for launching the solver
- `meshes` folder contains meshes that are needed by the solver
- `inputs` folder contains input parameters for setting up a test case for the Poisson_dg solver

Let us see how to set up a test cases in 2D. The file `inputs/circle.2D.arc` is used here for demonstration

#### Mesh ####

The mesh `circle_cut.msh` is provided in the `inputs/circle.2D.arc` file

```xml
  <meshes>
    <mesh>
      <filename>meshes/circle_cut.msh</filename>
    </mesh>
  </meshes>
```
Please not that use version 4.1 `.msh` file from `Gmsh`.

#### properties ###

The value of constant source term $\mathcal{f}$  can be provided in  `circle.2D.arc` file

```xml
  <fem>
    <f>5.5</f>
  </fem>
```


#### Boundary conditions ####

The Dirichlet boundary conditions  are provided in `circle.2D.arc` file

```xml
    <boundary-conditions>
      <dirichlet>
        <surface>horizontal</surface>
        <value>0.5</value>
      </dirichlet>
    </boundary-conditions>
```

So in the snippet above, three Dirichlet condition $u=0.5$ is  applied to border ('horizontal') which is a group of edges in the mesh file `circle_cut.msh`.

If needed, the Neumann  boundary conditions  can also be provided in `circle.2D.arc` file

```xml
    <boundary-conditions>
      <neumann>
        <surface>other</surface>
        <value>0.0</value>
      </neumann>
    </boundary-conditions>
```

#### Post Process ####

<img align="left" width="200" src="https://github.com/user-attachments/assets/66b9449e-e2f7-4607-b910-231def7d2f67" alt="poisson_1_large" />
<img align="left" width="200" src="https://github.com/user-attachments/assets/a68dd3d8-3f9f-424e-8d6a-33e34e41c04b" alt="poisson_1_large" />
Post processing is controlled via

```xml
  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>
```

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).
