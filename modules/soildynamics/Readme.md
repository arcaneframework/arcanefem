# Soildynamics
![alps](https://github.com/user-attachments/assets/878434f7-080c-4c6d-a470-b674df4b6996)

Here we deal with linear solid-mechanics governed by a system of PDE modeling the deformation of elastic bodies. The solver, here is a 2D unstructured mesh linear elasticity solver for soildynamics, which uses FEM to search for vector solution of displacement unknown $\mathbf{u}=(u_1,u_2)$, since transient $\mathbf{u}$ changes with time.

## Mathematics ##

#### Problem description ####

Under small elastic deformation, the steady  2D elastodynamics/soildynamics equation (linear elastic system) on  domain $\Omega$ reads

$$-\nabla\cdot\sigma(\mathbf{x}) +  \rho \mathbf{f}(\mathbf{x})=\rho \ddot{\mathbf{u}}(\mathbf{x}) \quad \mathbf{x}\in\Omega $$

here, $\sigma(\mathbf{x})$ is stress tensor, $\mathbf{f}(\mathbf{x})$ is the body force per unit volume, $\rho$ is the density of the material, and $\ddot{\mathbf{u}}(\mathbf{x})=\partial^2 \mathbf{u}(\mathbf{x})/\partial t^2$ is the acceleration for the displacement field $\mathbf{u}(\mathbf{x})$. The stress tensor  $\sigma(\mathbf{x})$ under isotropic elastic conditions is given by

$$ \sigma(\mathbf{x}) = \lambda(\nabla\cdot\mathbf{u}(\mathbf{x}))\mathbb{I} + \mu (\nabla\mathbf{u}(\mathbf{x}) + \left(\nabla\mathbf{u}(\mathbf{x})\right)^\text{T}) $$

here, $\lambda\in\mathbb{R}^{+}$ and $\mu\in\mathbb{R}^{+}$ are the Lame's elasticity parameters for the homogeneous material, $\mathbb{I}$ is the identity tensor, and $\mathbf{u}(\mathbf{x})$ is the displacement field vector. This governing PDE is also knows as Navier's equation.

#### Variational formulation ####

Without entering into the details, using test function $\mathbf{v}\in\mathbb{V}$ and trial function $\mathbf{u}\in\mathbb{V}$ with $\mathbb{V}$ being the suitable FE functional space,  the varaiational formulation for the soildynamics is more or less the same as that of the ealstodynamics with additional terms for boundary to absorb the outgoing wave, we do so using paraxial elements, the variational  problem reads:

$$\int_{\Omega} \rho \ddot{\mathbf{u}}(\mathbf{x})\cdot \mathbf{v}(\mathbf{x}) +  \int_{\Omega} \lambda \nabla \cdot \mathbf{u}(\mathbf{x}) \nabla \cdot \mathbf{v}(\mathbf{x}) + 2\mu\varepsilon(\mathbf{u}(\mathbf{x})):\varepsilon(\mathbf{u}(\mathbf{x})) +  \int_{\partial\Omega_P} \mathcal{P}(\mathbf{u}(\mathbf{x})) \cdot \mathbf{v}(\mathbf{x}) - \int_{\Omega}\rho\mathbf{f}(\mathbf{x})\cdot{\mathbf{v}(\mathbf{x})} - \int_{\partial\Omega_N} \mathbf{t}(\mathbf{x}) \cdot \mathbf{v}(\mathbf{x}) = 0 $$

here, $\mathbf{t}(\mathbf{x})$ is the traction vector imposed on Neumann boundary $\Omega_N$, $\mathcal{P}(\mathbf{u}(\mathbf{x}))$ is Paraxial term which is split into two terms (for the LHS and RHS) upon time discretization,  and  $\varepsilon(\mathbf{u}(\mathbf{x})) = \varepsilon_{ij}(\mathbf{u}(\mathbf{x}))$ is the strain tensor given by

$$\varepsilon_{ij}(\mathbf{u}) = \frac{1}{2}(\frac{\partial{u}_i}{\partial{x}_j} + \frac{\partial{u}_j}{\partial{x}_i} )$$

## The code ##

This XML configuration file is used for setting up an Soildynamics problem simulation in ArcaneFEM. Below is a detailed explanation of each section in the configuration for one such file `3d.double-couple.paraxial.soil.arc`.

###### Mesh Configuration

The mesh configuration section specifies the mesh file to be used in the simulation:

```xml
<meshes>
  <mesh>
    <filename>meshes/cube_double_couple_3d.msh</filename>
  </mesh>
</meshes>
```

- **Mesh File:** Defines the mesh file (`cube_double_couple_3d.msh`) to be used in the simulation. Note that this file should be compatible with version 4.1 `.msh` format from `Gmsh`.

###### FEM Configuration

The Finite Element Method (FEM) configuration is provided in the `Test.bar.arc`.

```xml
  <fem>
    <tmax>.1</tmax>
    <dt>0.01</dt>
    <cs>21.5</cs>
    <cp>40.</cp>
    <rho>9.0</rho>
    <paraxial-boundary-condition>
      <surface>leftsur</surface>
    </paraxial-boundary-condition>
    <double-couple>
      <north-node-name>dctop</north-node-name>
      <south-node-name>dcbot</south-node-name>
      <east-node-name>dcleft</east-node-name>
      <west-node-name>dcright</west-node-name>
      <method>force-based</method>
      <double-couple-input-file>data/force_loading_dc.txt</double-couple-input-file>
    </double-couple>
    <linear-system>
      <solver-backend>hypre</solver-backend>
    </linear-system>
  </fem>
```

Let us explain this point wise 

- **Time parameters:** The Maximum Time (tmax) time is set to `.1`. Time Step (dt) for the simulation is set to `0.01`. Time Discretization is set via Newmark-beta by default here.

  ```xml
  <tmax>.1</tmax>
  <dt>0.01</dt>
  ```

- **Material Properties:** The Density (rho) of the material is set to `9.0`. The primary wave speed `cp` is set to `40.` and the secondary speed `cs` is set to `20.1`..

  ```xml
    <cs>21.5</cs>
    <cp>40.</cp>
    <rho>9.0</rho>
  ```

- **Double Couple Sorce Condition:** Double couple is a way to introduce the earthquake source for the simulation. Within the mesh we would need to have four mesh nodes (groups) to which the double couple source is applied. The four points which are orthogonal should be provided appropiately, in our case the four points are named `dctop`, `dcbot`, `dcleft`, and `dcright`. In order to provide the transient loading to be applied at the double couple points, a text file `force_loading_dc.txt` is provided.  

  ```xml
    <double-couple>
      <north-node-name>dctop</north-node-name>
      <south-node-name>dcbot</south-node-name>
      <east-node-name>dcleft</east-node-name>
      <west-node-name>dcright</west-node-name>
      <method>force-based</method>
      <double-couple-input-file>data/force_loading_dc.txt</double-couple-input-file>
    </double-couple>
  ```

- **Paraxial Boundary Condition:** Defines the Paraxial boundary or in other words absorbing boundary condition for the seismic waves to exit the domain. Here we provide the mesh tag of the surface to which such conditions should be applied, `leftsur` in the case below.

  ```xml
    <paraxial-boundary-condition>
      <surface>leftsur</surface>
    </paraxial-boundary-condition>
  ```

- **Linear System Configuration:** Specifies the linear system settings, including the solver backend (`hypre`) which will use preconditioned Conjugate Gradient with AMG preconditioner.

  ```xml
  <linear-system>
    <solver-backend>hypre</solver-backend>
  </linear-system>
  ```


###### Post-Processing Configuration

The post-processing configuration is specified to control how and when results are saved:

```xml
  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
     <variable>V</variable>
     <variable>A</variable>
   </output>
  </arcane-post-processing>
```

- **Output Period:** Sets the interval at which results are saved.

- **Output Variables:** Lists the variables (`U` , `V`, and `A` ) which is the displacement, velocity, and acceleration vectors to be included in the output.



#### Post Process ####

For post processing the `ensight.case` file is outputted (in `output/depouillement` folder), which can be read by PARAVIS/ParaView. The output is of the $\mathbb{P}_1$ FE order (on nodes).

