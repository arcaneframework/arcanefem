# Aerodynamics with ArcaneFEM #
Potential Flow Theory Applied to Airfoil and Delta Wing Analysis is tacked in this module.

<img width="400" align="left" src="https://github.com/user-attachments/assets/6325a135-622a-4250-9d37-303898f35671" alt="Test_1_large_psi_new" />
<img width="400" align="left" src="https://github.com/user-attachments/assets/8c6cdf67-4cd1-4ebb-9bab-b556a0089fe0" alt="Test_1_large_psi_new" />
<img width="400" align="left" src="https://github.com/user-attachments/assets/03e3a781-cc8d-4d6a-933b-cf57ec746b5d" alt="Test_1_large_psi_new" />

Potential flow theory is a mathematical framework used to describe fluid flow under the assumptions that the fluid is inviscid, incompressible, and irrotational. These simplifications make the equations of motion more tractable, allowing for analytical solutions that provide insight into flow behavior.

In potential flow theory, the flow is characterized by a scalar function known as the **velocity potential**, $\phi$, which satisfies **Laplace’s equation**:

$$ \nabla^2 \phi = 0 $$

where $\nabla^2$ is the Laplacian operator. The velocity field $\mathbf{V}$ is derived as the gradient of the velocity potential:

$$ \mathbf{V} = \nabla \phi $$

Since the flow is irrotational, the **vorticity** $\mathbf{\omega}$ is zero:

$$ \mathbf{\omega} = \nabla \times \mathbf{V} = \nabla \times \nabla \phi = 0 $$

which ensures that the potential function $\phi$ is well-defined. In the case of **incompressible flow**, the continuity equation reduces to Laplace’s equation, reinforcing the applicability of potential flow analysis.

### Application to 2D and 3D Wing Profiles
Potential flow theory is applied to aerodynamic shapes such as the **NACA 0012 airfoil** and **3D delta wings** to determine flow characteristics like velocity distribution and pressure coefficients.

### Limitations and Practical Considerations
While potential flow theory provides a foundational understanding of aerodynamics, it neglects the effects of viscosity and turbulence, which are critical in real-world applications. Boundary layer separation, flow transition, and vortex shedding are phenomena not captured by potential flow theory but are essential for accurate aerodynamic predictions. Hence, potential flow analysis is often supplemented with computational fluid dynamics (CFD) or experimental testing for a more comprehensive evaluation.

Despite these limitations, potential flow remains a tool in preliminary aircraft design, lift estimation, and flow visualization.

## 3D simulation result of a UAV flow/pressure analysis ##
<img src="https://github.com/user-attachments/assets/8fb3bb53-6c7e-4fe9-b240-c97e1209a8a9" alt="Test_1_large_psi_new_new" style="zoom: 50%;" />
<img src="https://github.com/user-attachments/assets/2393ea17-7741-44d3-be4d-57ca9ec2b4ae" alt="Test_1_large_psi_new_new" style="zoom: 50%;" />

## 2D simulation result of a airfoil flow/pressure analysis ##
<img width="700" align="left" src="https://github.com/arcaneframework/arcanefem/assets/52162083/2c21cab5-5d7f-4bd9-a364-2b1f54e70edf" alt="Test_1_large_psi_new" />
<img src="https://github.com/arcaneframework/arcanefem/assets/52162083/8c691cee-d8e8-463a-b9b1-c00d016386f5" alt="Test_1_large_psi_new_new" style="zoom: 50%;" />
