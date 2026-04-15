# Solving a nonlinear Fourier equation with FEM and Arcane #

Here a nonlinear Fourier equation, that governs steady state heat conduction is solved using FEM in Arcane. The code here is a simple 2D unstructured mesh Galerkin FEM solver.

<img src="https://github.com/arcaneframework/arcanefem/assets/52162083/cf86f60f-360f-491b-a234-9631fc27af45" alt="Screenshot from 2022-12-19 16-25-59" style="zoom: 50%;" />


## Theory of nonlinear heat conduction ##

#### Problem description ####

The steady state 2D nonlinear heat conduction equation is solved for a closed meshed domain $\Omega^h$ in order to know the temperature $T(x,y)$ within the domain. The equation reads

$$\frac{\partial}{\partial x}\left(\lambda(T) \frac{\partial T}{\partial x} \right) + \frac{\partial}{\partial y}\left(\lambda(T) \frac{\partial T}{\partial y} \right)+ \dot{\mathcal{Q}} = 0  \quad \forall (x,y)\in\Omega^h $$

or in a more compact form

$$\nabla(\lambda(T) \nabla T) + \dot{\mathcal{Q}} = 0 \quad \forall (x,y)\in\Omega^h.$$

Here, $\lambda(T)$ is the temperature dependent thermal conductivity of the material and $\dot{\mathcal{Q}}$ is the heat generation source.
The nonlinearity in the above equation is due to the coefficient $\lambda(T)$. 

To complete the problem description, two first type (Dirichlet) boundary conditions are applied to this problem:

- $T = 0.0 &deg C \quad \forall(x,y)\in\partial\Omega^h_{\text{Gauche}}\subset\partial \Omega^h,$ and

- $T = 1.0 &deg C \quad \forall(x,y)\in\partial\Omega^h_{\text{Droite}}\subset\partial \Omega^h,$ 

in addition, all other boundaries $\partial \Omega^h_{N} = \partial \Omega^h \setminus (\partial\Omega^h_{\text{Gauche}} \cup \partial\Omega^h_{\text{Droite}})$ are exposed to second type (Neumann) boundary condition:
- $\lambda(T) \nabla T \cdot \mathbf{n}|_{\partial \Omega^h_{N} } = \overline{q} \cdot \mathbf{n}|_{\partial \Omega^h_{N} } = 0$

Finally, the heat-source term is set to zero

$\dot{\mathcal{Q}}=0$



### Finite element description of Fourier equation ###


In this case, the variational formulation in $H^1_{0}(\Omega) \subset H^1{\Omega}$ reads

search FEM trial function $u^h(x,y)$ satisfying

$$- \int_{\Omega^h}\lambda(u^h) \nabla u^h \nabla  v^h + \int_{\partial\Omega_N} (\overline{q} \cdot \mathbf{n}) v^h + \int_{\Omega^h}\dot{\mathcal{Q}} v^h = 0 \quad \forall v^h\in H^1_0(\Omega^h)$$

given

$u^h=0.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{Gauche}}$,

$u^h=1.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{Droite}}$ ,

$\int_{\Omega^h_{\text{Droite}}}(\mathbf{q} \cdot \mathbf{n}) v^h=15$,
    
$\int_{\Omega^h_{\text{Gauche}}}(\mathbf{q} \cdot \mathbf{n}) v^h=0$,

$\int_{\Omega^h}\dot{\mathcal{Q}} v^h=1\times10^5$, and


Please note that the above equation is often called as the weak formulation of the Fourier equation and in fact the finite element variable $u^h$ is an appoximation of temperature $T$.

## Problem Validation ##

### Thermal Conductivity ###
The inhomogenous thermal conductivity $\lambda(T)$ is a function of temperature, is given by

$$\lambda(T) = (1 + T)^m ,$$

such that setting $m=0$ brings us back to the homogenous $\lambda$ and linear Fourier problem.

###  Exact Solution ###
This definition along with the above BCs permit to obtain an analytical solution the nonlinear Fourier problem in Cartesian coordinates given by

$$T(x, y) = ((2^{m+1} - 1)x + 1)^{1/(1+m)} - 1 .$$

## The code ##

#### Heat Source ###

The value of heat source $\dot{\mathcal{Q}}$ can be provided in  `Test.nonlinear.conduction.arc` file

```xml
  <Fem1>
    <qdot>0.0</qdot>
  </Fem1>
```

#### Mesh ####

The mesh `unit_square.msh` is provided in the `Test.nonlinear.conduction.arc` file

```xml
  <meshes>
    <mesh>
      <filename>unit_square.msh</filename>
    </mesh>
  </meshes>
```

Please not that use version 4.1 `.msh` file from `Gmsh`.

#### Boundary conditions ####

The Dirichlet (constant temperature) boundary conditions  are provided in `Test.conduction.arc` file

```xml
    <dirichlet-boundary-condition>
      <surface>Cercle</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Bas</surface>
      <value>5.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Haut</surface>
      <value>21.0</value>
    </dirichlet-boundary-condition>
```

So in the snippet above, three Dirichlet conditions are applied ($50 \degree C, 5.0 \degree C, 21.0 \degree C$)  on three borders ('cercle', 'Bas', 'Haut').

The Neumann  boundary conditions  are also provided in `Test.conduction.arc` file

```xml
    <neumann-boundary-condition>
      <surface>Droite</surface>
      <value>15.0</value>
    </neumann-boundary-condition>
    <neumann-boundary-condition>
      <surface>Gauche</surface>
      <value>0.0</value>
    </neumann-boundary-condition>
```



#### Post Process ####

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).

#### multi-material example ####

<img src="https://github.com/arcaneframework/arcanefem/assets/52162083/eeac62de-3b5f-4264-a643-c6652a5693e8" alt="test_2_small_new"  />

#### manufactured solution example ####

![manu](https://github.com/arcaneframework/arcanefem/assets/52162083/5226528f-a5e1-4826-9b0b-36dc52ff57b3)


![manu](https://github.com/arcaneframework/arcanefem/assets/52162083/9237d686-2791-4852-b929-4d0c7e5f8df7)

