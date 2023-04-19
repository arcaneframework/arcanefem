# Transient elastodynamics

Here we deal with linear solid-mechanics governed by a system of PDE modeling the deformation of elastic bodies. The solver, here is a 2D unstructured mesh linear elasticity solver, which uses FEM to search for vector solution of displacement unknown $\mathbf{u}=(u_1,u_2)$, since transient $\mathbf{u}$ changes with time.



## Mathematics ##

#### Problem description ####

Under small elastic deformation, the steady  2D elastodynamics equation (linear elastic system) on  domain $\Omega$ reads

$$-\nabla\cdot\sigma(\mathbf{x}) +  \rho \mathbf{f}(\mathbf{x})=\rho \ddot{\mathbf{u}}(\mathbf{x}) \quad \mathbf{x}\in\Omega $$

here, $\sigma(\mathbf{x})$ is stress tensor, $\mathbf{f}(\mathbf{x})$ is the body force per unit volume, $\rho$ is the density of the material, and $\ddot{\mathbf{u}}(\mathbf{x})=\partial^2 \mathbf{u}(\mathbf{x})/\partial t^2$ is the acceleration for the displacement field $\mathbf{u}(\mathbf{x})$. The stress tensor  $\sigma(\mathbf{x})$ under isotropic elastic conditions is given by

$$ \sigma(\mathbf{x}) = \lambda(\nabla\cdot\mathbf{u}(\mathbf{x}))\mathbb{I} + \mu (\nabla\mathbf{u}(\mathbf{x}) + \left(\nabla\mathbf{u}(\mathbf{x})\right)^\text{T}) $$

here, $\lambda\in\mathbb{R}^{+}$ and $\mu\in\mathbb{R}^{+}$ are the Lame's elasticity parameters for the homogeneous material, $\mathbb{I}$ is the identity tensor, and $\mathbf{u}(\mathbf{x})$ is the displacement field vector. This governing PDE iis also knows as Navier's equation. 

#### Variational formulation ####

Without entering into the details, using test function $\mathbf{v}\in\mathbb{V}$ and trial function $\mathbf{u}\in\mathbb{V}$ with $\mathbb{V}$ being the suitable FE functional space,  the varaiational formulation for the elastodynamic problem reads:

$$\int_{\Omega} \rho \ddot{\mathbf{u}}(\mathbf{x})\cdot \mathbf{v}(\mathbf{x}) +  \int_{\Omega} \lambda \nabla \cdot \mathbf{u}(\mathbf{x}) \nabla \cdot \mathbf{v}(\mathbf{x}) + 2\mu\varepsilon(\mathbf{u}(\mathbf{x})):\varepsilon(\mathbf{u}(\mathbf{x})) - \int_{\Omega}\rho\mathbf{f}(\mathbf{x})\cdot{\mathbf{v}(\mathbf{x})} - \int_{\partial\Omega_N} \mathbf{t}(\mathbf{x}) \cdot \mathbf{v}(\mathbf{x}) = 0 $$

here, $\mathbf{t}(\mathbf{x})$ is the traction vector imposed on Neumann boundary $\Omega_N$, and  $\varepsilon(\mathbf{u}(\mathbf{x})) = \varepsilon_{ij}(\mathbf{u}(\mathbf{x}))$ is the strain tensor given by

$$\varepsilon_{ij}(\mathbf{u}) = \frac{1}{2}(\frac{\partial{u}_i}{\partial{x}_j} + \frac{\partial{u}_j}{\partial{x}_i} )$$  

## The code ##



#### Post Process ####

For post processing the `ensight.case` file is outputted, which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).
