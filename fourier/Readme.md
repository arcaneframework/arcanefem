# Solving Fourier equation with FEM and Arcane #

Here Fourier equation, that governs steady state heat conduction is solved using FEM in Arcane. The code here is a simple 2D unstructured mesh Galerkin FEM solver.

![Screenshot from 2022-12-19 16-25-59](https://github.com/arcaneframework/arcanefem/assets/52162083/cf86f60f-360f-491b-a234-9631fc27af45)


## Theory of heat conduction ##

#### Problem description ####

The steady state 2D heat conduction equation is solved for a closed meshed domain $\Omega^h$ in order to know the temperature $T(x,y)$ within the domain. The equation reads

$$\frac{\partial}{\partial x}\left(\lambda \frac{\partial T}{\partial x} \right) + \frac{\partial}{\partial y}\left(\lambda \frac{\partial T}{\partial y} \right)+ \dot{\mathcal{Q}} = 0  \quad \forall (x,y)\in\Omega^h $$

or in a more compact form

$$\nabla(\lambda\nabla T) + \dot{\mathcal{Q}} = 0 \quad \forall (x,y)\in\Omega^h.$$

Here, $\lambda$ is the thermal conductivity of the material and $\dot{\mathcal{Q}}$ is the heat generation source.



To complete the problem description,  three first type (Dirichlet) boundary conditions are applied to this problem:

$T = 50.0 \degree C \quad \forall(x,y)\in\partial\Omega^h_{\text{Cercle}}\subset\partial \Omega^h,$

$T = 5.0\degree C \quad \forall(x,y)\in\partial\Omega^h_{\text{Bas}}\subset\partial \Omega^h,$ and

$T  = 21.0\degree C  \quad \forall(x,y)\in\partial\Omega^h_{\text{Haut}}\subset\partial \Omega^h,$

in addition, other boundaries $\partial\Omega^h_N$ are exposed to  second type (Neumann) boundary condition:

- first Neumann condition $\partial\Omega^h_{\text{Droite}}$ the derivative of temperature (heat flux) is non-null - influx boundary

$$\mathbf{q}\cdot\mathbf{n}|_{\partial \Omega^h_{\text{Droite}}} = 15.0$$

- second Neumann condition $\partial\Omega^h_{\text{Gauche}}$ the derivative of temperature (heat flux) is null - insulation boundary

$$\mathbf{q}\cdot\mathbf{n}|_{\partial \Omega^h_{\text{Gauche}}} = 0$$

Finally a uniform heat-source is present within the domain

$\dot{\mathcal{Q}}=1\times10^5$



### Finite element description of Fourier equation ###



We work with approximation, $\lambda$ is homogeneous $\lambda : \Omega^h \in \mathbb{R}^{+}$, in this case  the variational formulation in $H^1_{0}(\Omega) \subset H^1{\Omega}$  reads

search FEM trial function $u^h(x,y)$ satisfying

$$- \int_{\Omega^h}\lambda\nabla u^h \nabla  v^h + \int_{\partial\Omega_N} (\overline{q} \cdot \mathbf{n}) v^h + \int_{\Omega^h}\dot{\mathcal{Q}} v^h = 0 \quad \forall v^h\in H^1_0(\Omega^h)$$

given

$u^h=50.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{Cercle}}$,

$u^h=5.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{Bas}}$ ,

$u^h=20.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{Haut}}$,

$\int_{\Omega^h_{\text{Droite}}}(\mathbf{q} \cdot \mathbf{n}) v^h=15$,

$\int_{\Omega^h_{\text{Gauche}}}(\mathbf{q} \cdot \mathbf{n}) v^h=0$,

$\int_{\Omega^h}\dot{\mathcal{Q}} v^h=1\times10^5$, and

$\lambda=1.75$

Please note that the above equation is often called as the weak formulation of the Fourier equation and in fact the finite element variable $u^h$ is an appoximation of temperature $T$. 

## The code ##

#### Thermal Conductivity ###

The value of thermal conductivity $\lambda$  and heat source $\dot{\mathcal{Q}}$ can be provided in  `Test.conduction.arc` file

```xml
  <Fem1>
    <lambda>1.75</lambda>
    <qdot>1e5</qdot>
  </Fem1>
```

#### Mesh ####

The mesh `plancher.msh` is provided in the `Test.conduction.arc` file

```xml
  <meshes>
    <mesh>
      <filename>plancher.msh</filename>
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

For post processing the `ensight.case` file is outputted, which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).

#### multi-material example ####

![test_2_large_new](https://github.com/arcaneframework/arcanefem/assets/52162083/f5a61339-3f9a-40c1-bc94-4b9090905e20)



