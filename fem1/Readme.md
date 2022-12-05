# Notes on the solver #

The code here is a simple FEM code used to solve a conduction problem on unstructured 2D mesh. 



## Theory ##

#### Problem description ####

The steady state 2D heat conduction equation is solved for a closed meshed domain $\Omega^h$ in order to know the temperature $T(x,y)$ within the domain. The equation reads

$$\frac{\partial}{\partial x}\left(\lambda \frac{\partial T}{\partial x} \right) + \frac{\partial}{\partial y}\left(\lambda \frac{\partial T}{\partial y} \right)+ \dot{\mathcal{Q}} = 0  \quad \forall (x,y)\in\Omega^h $$

or in a more compact form

$$\nabla(\lambda\nabla T) + \dot{\mathcal{Q}} = 0 \quad \forall (x,y)\in\Omega^h.$$

Here, $\lambda$ is the thermal conductivity of the material and $\dot{\mathcal{Q}}$ is the heat generation source.



To complete the problem description,  three first type (Dirichlet) boundary conditions are applied to this problem

$$T|_{\partial\Omega^h_{\text{Cercle}}} = 50.0 \degree C, \quad T|_{\partial\Omega^h_{\text{Bas}}} = 5.0\degree C, \quad \text{and}\quad T|_{\partial\Omega^h_{\text{Haut}}} = 21.0\degree C$$,

in addition, other boundaries $\partial\Omega^h_N$ are exposed to  second type (Neumann) boundary condition, the derivative of temperature (heat flux) is zero - insulation boundary

$$\mathbf{q}\cdot\mathbf{n}|_{\partial \Omega^h_N} = 0$$



We work with approximation, $\lambda$ is homogeneous $\lambda : \Omega^h \in \mathbb{R}^{+}$, there is not heat source present $\dot{\mathcal{Q}}=0$.  The variational formulation in $H^1_{0}(\Omega) \subset H^1{\Omega}$  reads

search FEM trial function $u^h(x,y)$ satisfying

$$- \int_{\Omega^h}\lambda\nabla u^h \nabla  v^h + \int_{\partial\Omega_N} (\overline{q} \cdot \mathbf{n}) v^h + \int_{\Omega^h}\dot{\mathcal{Q}} v^h = 0 \quad \forall v^h\in H^1_0(\Omega^h)$$

given

  $$u^h=50.0 \quad \forall (x,y)\in\Omega^h_{\text{Cercle}}, $$

 $u^h=5.0 \quad \forall (x,y)\in\Omega^h_{\text{Bas}}$ ,

$$u^h=20.0 \quad \forall (x,y)\in\Omega^h_{\text{Haut}}$$,

 $\int_{\Omega^h_N}(\mathbf{q} \cdot \mathbf{n}) v^h=0$,

$\int_{\Omega^h}\dot{\mathcal{Q}} v^h=0$, and

$\lambda=1.75$



## The code ##

#### Thermal Conductivity ###

The value of thermal conductivity $\lambda$ can be provided in  `FemTest1.arc` file

```xml
  <Fem1>
    <lambda>1.75</lambda>
  </Fem1>
```

#### Mesh #### 

The mesh `plancher.msh` is provided in the `FemTest1.arc` file 

```xml
  <meshes>
    <mesh>
      <filename>plancher.msh</filename>
    </mesh>
  </meshes>
```

Please not that use version 4 `.msh` file from `Gmsh`. 

#### Boundary conditions ####

These are setup via the `Fem1Module::_initBoundaryconditions()` within the `Fem1Module.cc`:

```c++
void Fem1Module::
_initBoundaryconditions()
{
  info() << "Init boundary conditions...";
    
  _applyOneBoundaryCondition("Cercle", 50.0);
  _applyOneBoundaryCondition("Bas", 5.0);
  _applyOneBoundaryCondition("Haut", 21.0);
}
```

the `_applyOneBoundaryCondition()` function is used to input the three Dirichlet boundary conditions here. For example, via  `_applyOneBoundaryCondition("Bas", 5.0);` we impose $5.0$ $\degree C$ on the border "Bas" (french for bottom). Note, "Bas" is a physical tag in the `plancher.msh` file.  So in the snippet above, three Dirichlet conditions are applied ($50 \degree C, 5.0 \degree C, 21.0 \degree C$)  on three borders ('cercle', 'Bas', 'Haut').

#### Zero flux and zero source condition ###

In FEM the source therm and the flux term are a part of RHS vector. Since here $\int_{\Omega^h}\dot{\mathcal Q} v^h = 0  $ and $\int_{\partial\Omega^h_N} ({\mathbf{q}} \cdot \mathbf{n}) v^h = 0$, the RHS contribution due to these terms is null. These are setup in `Fem1Module.cc`

 ```c++
 void Fem1Module::
 _computeGeneralizedFluxes()
 {
   // TODO: Loop over all faces on the border instead
   m_rhs_vector.fill(0.0);
 }
 ```

and

```c++
void Fem1Module::
_computeSourceTerm()
{
  // TODO: Loop over all cells and fill the source term
  m_rhs_vector.fill(0.0);
}
```



#### Post Process ####

For post processing the `ensight.case` file is outputted, which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).
