# Solving heat conduction equation with FEM and Arcane #

Here time dependent parabolic problem of heat conduction, is solved using FEM in Arcane.

## Theory of heat conduction ##

#### Problem description ####

The transient 2D heat conduction equation is solved for a closed meshed domain $\Omega^h$ in order to know the temperature $T(x,y)$ within the domain. The equation reads

$$\frac{\partial T}{\partial t} - \frac{\partial}{\partial x}\left(\lambda \frac{\partial T}{\partial x} \right) - \frac{\partial}{\partial y}\left(\lambda \frac{\partial T}{\partial y} \right)+ \dot{\mathcal{Q}} = 0  \quad \forall (x,y)\in\Omega^h $$

or in a more compact form

$$\partial_tT - \nabla(\lambda\nabla T) + \dot{\mathcal{Q}} = 0 \quad \forall (x,y)\in\Omega^h.$$

Here, $\lambda$ is the thermal conductivity of the material and $\dot{\mathcal{Q}}$ is the heat generation source.



To complete the problem description,  one first type (Dirichlet) boundary conditions are applied to this problem:

$T = 10.0 \degree C \quad \forall(x,y)\in\partial\Omega^h_{\text{left}}\subset\partial \Omega^h,$

in addition, we are provided with initial conditions:

$T(x,y,0)=30.0$

Finally no  heat-source is present within the domain 

$\dot{\mathcal{Q}}=0$



#### FEM discretization

We work with approximation, $\lambda$ is homogeneous $\lambda : \Omega^h \in \mathbb{R}^{+}$, and apply the implicit Euler finite difference approximation fir time,  in this case  the variational formulation in $H^1_{0}(\Omega) \subset H^1{\Omega}$  reads

search FEM trial function $u_n^h(x,y)$ satisfying

$$- \int_{\Omega^h} \frac{u^h_n - u^h_{n-1}}{\delta t}+ \int_{\Omega^h}\lambda\nabla u_n^h \nabla  v^h  = 0 \quad \forall v^h\in H^1_0(\Omega^h)$$

given

$u_n^h=10.0 \quad \forall (x,y)\in\partial\Omega^h_{\text{left}}$,

$u_0^h=30.0 \quad \forall (x,y)\in\Omega^h$, and

$\lambda=1.75$



We also are provided with  $t_{max}=40$ seconds and $\delta t=0.4$, as such the simulation is of 100 time-steps.



## The code ##

#### Thermal Conductivity ###

The value of thermal conductivity $\lambda$  and heat source $t_{max}, \delta t$, and $T_{init}$ can be provided in  `Test.conduction.arc` file

```xml
  <Fem1>
    <lambda>1.75</lambda>
    <tmax>40.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
  </Fem1>
```

#### Mesh #### 

The mesh `plate.msh` is provided in the `Test.conduction.arc` file 

```xml
  <meshes>
    <mesh>
      <filename>plate.msh</filename>
    </mesh>
  </meshes>
```

Please not that use version 4.1 `.msh` file from `Gmsh`. 

#### Boundary conditions ####

The Dirichlet (constant temperature) boundary conditions  are provided in `Test.conduction.arc` file

```xml
    <dirichlet-boundary-condition>
      <surface>Left</surface>
      <value>10.0</value>
    </dirichlet-boundary-condition>
```

So in the snippet above, three Dirichlet conditions are applied $10 \degree C$  on  border tagged as 'left' in the .



#### Post Process ####

For post processing the `ensight.case` file is outputted, which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes). In the `Test.conduction.arc` file we provide

```xml
  <arcane-post-processing>
   <output-period>2</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>
```

This means we would like to postprocess`NodeTemprature` variable which is the FEM variable $u^h_n$ and physical variable $T$.  Furthermore the `<output-period>2</output-period>` indicates that for every 2 time iteration dump simulation results in the post-processing file. So for this simulation which has 100 time steps we will post-process 50 times. 
