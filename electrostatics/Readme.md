# Electrostatics  with ArcaneFEM #

The code here is a simple 2D unstructured mesh Galerkin FEM solver for Electrostatics.

![electrostatics_2](https://github.com/arcaneframework/arcanefem/assets/52162083/959988a3-1717-4449-b412-14cbd1582367)

## Theory ##

#### Problem description ####

Assuming no current and steady-state change distribution, 2D Electrostatic equation is solved for a closed meshed domain $\Omega^h$ in order to know the electric field $\mathbf{E}(x,y)$ within the domain. The equation reads

$$\text{div}\mathbf{E} = \rho/\epsilon   \quad \forall (x,y)\in\Omega^h $$

$$\text{curl}\mathbf{E} = 0   \quad \forall (x,y)\in\Omega^h $$

here, $\rho$ and $\epsilon$ are charge density and permittivity of free space, respectively.

Introducing, electrostatic potential $\phi$ such that

$$\mathbf{E}=-\nabla\phi$$

we get the following form of equation

$$-\nabla\phi=\rho/\epsilon  \quad \forall (x,y)\in\Omega^h$$

To complete the problem description,   first type (Dirichlet) boundary conditions is applied to this problem, we assume that external border is held at zero potential:

$\phi = 0.0 \quad \forall(x,y)\in\partial\Omega^h_{\text{external}}\subset\partial \Omega^h,$

and internally there are two rods held at +1 and -1 volt receptively,  which is again translated to two Dirichlet boundary conditions

$\phi = 1.0 \quad \forall(x,y)\in\partial\Omega^h_{\text{rod1}}\subset\partial \Omega^h,$

$\phi = -1.0 \quad \forall(x,y)\in\partial\Omega^h_{\text{rod2}}\subset\partial \Omega^h,$

### Finite Element method ###

In this case  the FEM variational formulation in $H^1_{0}(\Omega) \subset H^1{\Omega}$  reads

search FEM trial function $\phi^h(x,y)$ satisfying

$$- \int_{\Omega^h}\nabla \phi^h \nabla  v^h + \int_{\partial\Omega_N} (\overline{q} \cdot \mathbf{n}) v^h + \int_{\Omega^h}(\rho/\epsilon) v^h = 0 \quad \forall v^h\in H^1_0(\Omega^h)$$

given

$\phi^h=0.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{external}}$,

$\phi^h=1.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{rod1}}$,

$\phi^h=-1. \quad \forall (x,y)\in\partial\Omega^h_{\text{rod2}}$,

$\int_{\Omega^h_{\text{N}}}(\mathbf{q} \cdot \mathbf{n}) v^h=0$ since no Neumann BC is present,

## The code ##

#### properties ###

The value of constant material property terms $\rho,\epsilon$  can be provided in  `Test.Electrostatics.arc` file

```xml
  <fem>
    <rho>0</rho>
    <epsilon>1</epsilon>
  </fem>
```

#### Mesh ####

The mesh `box-rods.msh` is provided in the `Test.Electrostatics.arc` file

```xml
  <meshes>
    <mesh>
      <filename>box-rods.msh</filename>
    </mesh>
  </meshes>
```

Please not that use version 4.1 `.msh` file from `Gmsh`.

#### Boundary conditions ####

The Dirichlet boundary conditions  are provided in `Test.Electrostatics.arc` file

```xml
    <dirichlet-boundary-condition>
      <surface>external</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>rod1</surface>
      <value>1.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>rod2</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
```

So in the snippet above, three Dirichlet condition $\phi=0$ is  applied to border ('external') which is a group of edges in the mesh file `box-rods.msh`. Similar logic applies to other two Dirichlet conditions.

If needed, the Neumann  boundary conditions  can also be provided in `Test.Electrostatics.arc` file

```xml
    <neumann-boundary-condition>
      <surface>Neumann</surface>
      <value>0.0</value>
    </neumann-boundary-condition>
```



#### Post Process ####

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).


## Gallery ##
- E-type interdigital capacitor
![image](https://github.com/arcaneframework/arcanefem/assets/52162083/822a3b8c-5a55-4a5b-8450-74224d2a257d)

![image](https://github.com/arcaneframework/arcanefem/assets/52162083/f38a3390-56ad-498d-867d-d830116957e3)
