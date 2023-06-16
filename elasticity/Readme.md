# linear elasticity
<img width="500" align="left" src="https://github.com/arcaneframework/arcanefem/assets/52162083/eb970ece-5fd3-4862-9b93-e8930a103ae9" />


Here, we focus on linear solid mechanics, which involves analyzing the behavior of elastic bodies through a system of partial differential equations (PDEs) that govern their deformation. The solver at hand is specifically designed for 2D unstructured mesh problems and operates as a linear elasticity solver utilizing the finite element method (FEM).

The primary objective of the solver is to determine the vector solution for the unknown displacements within the elastic body. By utilizing FEM, it efficiently searches for a solution that satisfies the given boundary conditions and accurately captures the deformation characteristics of the system. This solver serves as a valuable tool for investigating and understanding the mechanical behavior of elastic materials in various engineering applications.


## Mathematics ##

#### Problem description ####

Under small elastic deformation, the steady  2D elastic deformation equation (linear elastic system) on  domain $\Omega$ reads

$$-\nabla\cdot\sigma(\mathbf{x})=\mathbf{f}(\mathbf{x}) \quad \mathbf{x}\in\Omega $$

here, $\sigma(\mathbf{x})$ is stress tensor, $\mathbf{f}(\mathbf{x})$ is the body force per unit volume. The stress tensor  $\sigma(\mathbf{x})$ under isotropic elastic conditions is given by

$$ \sigma(\mathbf{x}) = \lambda(\nabla\cdot\mathbf{u}(\mathbf{x}))\mathbb{I} + \mu (\nabla\mathbf{u}(\mathbf{x}) + \left(\nabla\mathbf{u}(\mathbf{x})\right)^\text{T}) $$

here, $\lambda\in\mathbb{R}^{+}$ and $\mu\in\mathbb{R}^{+}$ are the Lame's elasticity parameters for the homogeneous material, $\mathbb{I}$ is the identity tensor, and $\mathbf{u}(\mathbf{x})$ is the displacement field vector. This governing PDE iis also knows as Navier's equation.

#### Variational formulation ####

Without entering into the details, the variational formulation for the Navier's equation reads

$$\int_{\Omega} \lambda \nabla \cdot \mathbf{u}(\mathbf{x}) \nabla \cdot \mathbf{v}(\mathbf{x}) + 2\mu\varepsilon(\mathbf{u}(\mathbf{x})):\varepsilon(\mathbf{u}(\mathbf{x})) - \int_{\Omega}\mathbf{f}(\mathbf{x})\cdot{\mathbf{v}(\mathbf{x})} - \int_{\partial\Omega_N} \mathbf{t}(\mathbf{x}) \cdot \mathbf{v}(\mathbf{x}) = 0 $$

here, $\mathbf{t}(\mathbf{x})$ is the traction vector imposed on Neumann boundary $\Omega_N$, and  $\varepsilon(\mathbf{u}(\mathbf{x})) = \varepsilon_{ij}(\mathbf{u}(\mathbf{x}))$ is the strain tensor given by

$$\varepsilon_{ij}(\mathbf{u}) = \frac{1}{2}(\frac{\partial{u}_i}{\partial{x}_j} + \frac{\partial{u}_j}{\partial{x}_i} )$$

## The code ##



#### Post Process ####

For post processing the `ensight.case` file is output, which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).
