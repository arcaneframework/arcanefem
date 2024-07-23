# Bilaplacian with ArcaneFEM #


<img width="500" align="left" src="https://github.com/arcaneframework/arcanefem/assets/52162083/9f183f44-cc7c-40cb-9b6b-8fefdf0f94bf"/>


The code here is a simple FEM code used to solve a bilaplacian problem on unstructured 2D mesh. Here bilaplacian refers to solution of two Laplace equations (hence the name bilaplacian) coupled through a Poisson equation.



## Theory ##

#### Problem description ####

The steady state 2D bilaplacian equation is solved for a closed square meshed domain $\Omega^h = (0,1)^2$ in order to know the vectorial unknowns $\{u_i(x,y)\}_{i=1}^2$, i.e, $u_1$ and $u_2$ within the domain. The system of equations read

$$\triangle u_1 + u_2  = 0  \quad \forall (x,y)\in\Omega^h $$

$$\triangle u_2  = f  \quad \forall (x,y)\in\Omega^h $$

There are no Neumann boundrary conditions attached to this problem, however to complete the problem description,  four first type (Dirichlet) boundary conditions are applied to this problem:

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Top}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Left}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Right}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Bot}}\subset\partial \Omega^h,$


This system of equations arises in various applications, such as fluid mechanics, electromagnetics, and elasticity. In fluid mechanics, for instance, it can be used to model the flow of an incompressible fluid over a flat plate. In this case, $u_1$ represents the stream function, which describes the flow velocity in the $x$-direction, while $u_2$ represents the velocity potential, which describes the flow velocity in the $y$-direction.

#### Variational formulation


In this case  the variational formulation we use the subset of square integrable Sobolev functional space   $H^1_{0}(\Omega) \subset H^1{\Omega}$. The FEM formulation then reads:

search vectorial FEM trial function $(u^h_1,u^h_2)\in\left[H^1_0(\Omega^h)\right]^2$ satisfying

$$ \int_{\Omega^h}\nabla u^h_1 \nabla  v_2^h +  \int_{\Omega^h}\nabla u^h_2 \nabla  v_1^h + \int_{\Omega^h} u^h_2   v_2^h + \int_{\Omega^h}f v_1^h = 0 \quad \forall (v_1^h,v_2^h)\in \left[H^1_0(\Omega^h)\right]^2$$

given:

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Top}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Left}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Right}}\subset\partial \Omega^h,$

$u_1 = 0.0  \quad \forall(x,y)\in\partial\Omega^h_{\text{Bot}}\subset\partial \Omega^h,$

$\int_{\Omega^h}f v_1^h=1\times10^5$



## The code ##

This XML configuration file (e.g `Test.direct.arc`) is used for setting up a Finite Element Method (FEM) simulation in ArcaneFEM. Below is a detailed breakdown of each section in the configuration.

##### Mesh Configuration ######

The mesh configuration section specifies the mesh file to be used in the simulation:

```xml
<meshes>
  <mesh>
    <filename>bilap.msh</filename>
  </mesh>
</meshes>
```

- Defines the mesh file (`bilap.msh`) to be used in the simulation. Note that this file should be compatible with version 4.1 `.msh` format from `Gmsh`.

###### FEM Configuration

```xml
<fem>
  <f>-1.0</f>
  <enforce-Dirichlet-method>WeakPenalty</enforce-Dirichlet-method>
  <penalty>1.e30</penalty>
  <dirichlet-boundary-condition>
    <surface>boundary</surface>
    <value>0.05</value>
  </dirichlet-boundary-condition>
  <linear-system name="SequentialBasicLinearSystem">
    <solver-method>direct</solver-method>
  </linear-system>
</fem>
```

- **Source Term (f):** The source term in the partial differential equation (PDE) is set to `-1.0`.
- **Dirichlet Method:** Specifies the method (`WeakPenalty`) for enforcing Dirichlet boundary conditions. And we set the penalty parameter (`1.e30`) for weak enforcement of boundary conditions.
- **Dirichlet Boundary Condition:** Defines the boundary condition on the specified surface (`boundary`) with a given value (`0.05`).
- **Linear System Configuration:** Specifies the linear system settings, including the solver method (`direct`) to be used for solving the FEM problem.



#### Post Process ####

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).

#### Tests available in this module ####

The tests are present in the form of `.arc` files with a prefix `Test.`:

| Name         | Dimension | Boundary Condition                     | Solver               | Comment                           |
| ------------ | --------- | -------------------------------------- | -------------------- | --------------------------------- |
| direct       | 2D        | Dirichlet only<br />Homogeneous Source | Sequential Direct LU | - Test Weak Penalty method for BC |
| internal_pcg | 2D        | Dirichlet only<br />Homogeneous Source | Arcane's PCG solver  |                                   |
|              |           |                                        |                      |                                   |

