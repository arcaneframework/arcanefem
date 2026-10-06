# Elastoplasticity

<img width="1592" height="364" alt="vp" src="https://github.com/user-attachments/assets/a91ffcd3-1ee2-4946-87c7-773c22c944d9" />

We focus on the nonlinear solid mechanics, specifically analyzing the elastoplastic response of solid bodies through a system of partial differential equations (PDEs) and internal state variable updates. The nonlinear solver implemented using Newton method resolves the  elastoplasticity problem utilizing the finite element method (FEM).

The primary objective of the solver is to determine the vector displacement field and path-dependent stress state within the body under applied loading. By utilizing FEM with an iterative Newton scheme and local return-mapping algorithms, it accurately captures non-reversible plastic deformation and material yielding. This solver serves as a valuable tool for investigating rate-independent plastic collapse and stress redistribution in engineering applications.

## Mathematics

#### 2D Problem description

We demosntrate the module with an example case that addresses the incremental analysis of a nonlinear elasto-plastic Von‑Mises material in 2D ($x$‑$y$). The domain of interest is a quarter‑cylinder with external radius $R_e$ and internal radius $R_i$.

Symmetry conditions are applied at the bottom ($y=0$) and left ($x=0$) boundaries. A uniform internal pressure $q$ is applied at the internal boundary. This pressure increases from 0 to

$$
q_{\text{lim}} = \frac{2}{\sqrt3} \sigma_0 \,\ln\!\left(\frac{R_e}{R_i}\right)
$$

which is the analytical collapse load for a perfectly plastic material (no hardening).


Under small-strain assumptions, the static equilibrium equation on the 2D domain $\Omega$ reads:

$$
-\nabla\cdot\sigma(\mathbf{x})=F_{\text{ext}}(\mathbf{x}) \quad \mathbf{x}\in\Omega
$$

where $\sigma(\mathbf{x})$ is the Cauchy stress tensor and $F_{\text{ext}}$ is the external body force per unit volume. The total strain tensor $\varepsilon(\mathbf{u})$ is split additively into elastic and plastic components:

$$
\varepsilon(\mathbf{u}) = \frac{1}{2}\left(\nabla\mathbf{u} + (\nabla\mathbf{u})^\text{T}\right) = \varepsilon^e + \varepsilon^p
$$

The stress tensor $\sigma(\mathbf{x})$ is related to the elastic strain via Hooke's law:

$$
\sigma(\mathbf{x}) = \lambda \operatorname{tr}(\varepsilon^e)\mathbb{I} + 2\mu \varepsilon^e
$$

where $\lambda\in\mathbb{R}^{+}$ and $\mu\in\mathbb{R}^{+}$ are the Lamé parameters for the homogeneous isotropic material, and $\mathbb{I}$ is the 2D identity tensor. Yielding is governed by the isotropic von Mises criterion:

$$
f(\sigma) = \sigma^{\text{eq}} - \sigma_0 - H p \le 0, \quad \sigma^{\text{eq}} = \sqrt{\frac{3}{2} \mathbf{s}:\mathbf{s}}
$$

where $\mathbf{s} = \sigma - \frac{1}{3}\operatorname{tr}(\sigma)\mathbb{I}$ is the deviatoric stress tensor, $\sigma_0$ is the initial yield strength, $H$ is the linear isotropic hardening modulus, and $p$ is the accumulated equivalent plastic strain.


#### Variational formulation and incremental solution

Without entering into the details, the variational formulation for the above description reads:

$$
\int_{\Omega} \sigma(\mathbf{u}(\mathbf{x})) : \varepsilon(\mathbf{v}(\mathbf{x})) \, \text{d}\Omega - \int_{\Omega} F_{\text{ext}}(\mathbf{x})\cdot\mathbf{v}(\mathbf{x}) \, \text{d}\Omega - \int_{\partial\Omega_N} \mathbf{t}(\mathbf{x}) \cdot \mathbf{v}(\mathbf{x}) \, \text{d}\Gamma = 0
$$

where $\mathbf{t}(\mathbf{x})$ is the surface traction vector applied on the Neumann boundary $\partial\Omega_N$.

Because the material response is path-dependent, the problem is solved incrementally over time/load steps $\Delta t$. Within step $n+1$, given state variables $(\sigma_n, \varepsilon^p_n, p_n)$ and a strain increment $\Delta\varepsilon = \varepsilon(\mathbf{u}_{n+1}) - \varepsilon(\mathbf{u}_n)$, the local stress update is computed using a radial return-mapping predictor-corrector strategy:

1. **Trial elastic state:**
   $$
   \sigma_{\text{elas}} = \sigma_n + \mathbf{C}^e : \Delta\varepsilon
   $$
   $$
   \mathbf{s}_{\text{elas}} = \sigma_{\text{elas}} - \frac{1}{3}\operatorname{tr}(\sigma_{\text{elas}})\mathbb{I}, \quad \sigma_{\text{elas}}^{\text{eq}} = \sqrt{\frac{3}{2} \mathbf{s}_{\text{elas}} : \mathbf{s}_{\text{elas}}}
   $$

2. **Yield condition check:**
   $$
   f_{\text{elas}} = \sigma_{\text{elas}}^{\text{eq}} - \sigma_0 - H p_n
   $$

3. **Plastic update:**
   If $f_{\text{elas}} \le 0$, the response is elastic ($\Delta p = 0$, $\sigma_{n+1} = \sigma_{\text{elas}}$).  
   Otherwise, plastic flow occurs and the increment in plastic strain is evaluated analytically as:
   $$
   \Delta p = \frac{f_{\text{elas}}}{3\mu + H}
   $$

4. **Stress & state variable correction:**
   $$
   \sigma_{n+1} = \sigma_{\text{elas}} - 3\mu \Delta p \frac{\mathbf{s}_{\text{elas}}}{\sigma_{\text{elas}}^{\text{eq}}}
   $$
   $$
   p_{n+1} = p_n + \Delta p, \quad \varepsilon^p_{n+1} = \varepsilon^p_n + \frac{3}{2} \Delta p \frac{\mathbf{s}_{\text{elas}}}{\sigma_{\text{elas}}^{\text{eq}}}
   $$

In the above procedure, steps 1 to 3 iteratively solved using a Newton method until the trial stress is consistent with the new state (indicated by convergence of Newton solve). The Newton algorithmic uses a elastoplastic tangent tensor (Jacobian in the Newton solve) 
$$
\mathbb{C}_{\text{tang}}^{\text{alg}} = \mathbb{C} - 3\mu \left( \frac{3\mu}{3\mu + H} - \beta \right) \mathbf{\mathit{n}}_{\text{elas}} \otimes \mathbf{\mathit{n}}_{\text{elas}} - 2\mu\beta \mathbb{D}\text{ev}.
$$
Consequently, the RHS now has an internal force 
$$
-F_{\text{int}} = - \int_{\Omega} \sigma_{n} : \varepsilon(\mathbf{v}(\mathbf{x}))
$$ 
term which a forms residual of the nonlinear. Therefore, a Newton step $k$ solves the linear system
$$
  \mathbb{C}_{\text{tang}}^{\text{alg}} \Delta u_k = F_{\text{ext}} - F_{\text{int}} + \texttt{traction\_BC}
$$ 
where $\Delta u_k = u_{k+1} - u_{k}$

#### NOTE: The  processes from 1 to 3 of getting $\mathbb{C}_{\text{tang}}^{\text{alg}}$ from the previous know states could be handled by ```MFront```


## Von Mises law through MGIS (MFrontGenericInterfaceSupport)

The constitutive law `VonMisesMGIS` delegates the integration of the von Mises
elastoplastic behaviour (with linear isotropic hardening) to the
[MGIS](https://thelfer.github.io/mgis/web/index.html) library, using a
behaviour generated by
[MFront](https://thelfer.github.io/tfel/web/index.html). The native
implementation of the law (`VonMises`) is kept and both laws can be compared
on the same case: they reproduce each other to machine precision.

#### Requirements

The `VonMisesMGIS` law is enabled at configure time when both the MGIS library
and the `mfront` code generator are found, for instance:

```bash
cmake <arcanefem-sources> -DMGIS_DIR=<mgis-installation-prefix> \
      -DCMAKE_PROGRAM_PATH=<tfel-installation-prefix>/bin
```

A `FindMFrontGenericInterface.cmake` module is provided: it searches the MGIS
package configuration file (usually `<prefix>/share/mgis/cmake`) from the
`MGIS_DIR` or `MFrontGenericInterface_DIR` hints. If MGIS or `mfront` is not
found, the `VonMisesMGIS` law is disabled and selecting it leads to a fatal
error, while the rest of the module is unaffected.

#### The MFront behaviour

The file [mfront/VonMises.mfront](mfront/VonMises.mfront) describes the
behaviour with the `IsotropicPlasticMisesFlow` domain specific language: the
yield function is $f(\sigma^{eq},p) = \sigma^{eq} - s_0 - H p$ where $s_0$ is
the initial yield strength and $H$ the linear isotropic hardening modulus. The
Young modulus, the Poisson ratio, the yield strength and the hardening slope
are material properties provided at runtime by the solver. The behaviour
sources are generated by `mfront` at build time and compiled as the shared
library `libVonMises.so` next to the executable.

#### Usage

The law is selected in the arc file with the same material parameters as the
native law (the hardening slope `H` defaults to the value derived from `E` as
in the native law, i.e. $H = EE_t/(E-E_t)$ with $E_t = E/100$):

```xml
<constitutive-law>
  <law>VonMisesMGIS</law>
  <von-mises-mgis>
    <library>libVonMises.so</library>
    <behaviour>VonMises</behaviour>
    <E>70.0e3</E>
    <nu>0.3</nu>
    <sig0>250.</sig0>
  </von-mises-mgis>
</constitutive-law>
```

The `library` option is the path of an MFront behaviour library (a path
relative to the current directory is resolved explicitly) and `behaviour` the
name of the behaviour inside this library. Any MFront-generated small strain
von Mises behaviour with the same material properties can be used.

The integration is performed on the CPU at each Gauss point of the Tria3,
Quad4, Quad8 and Quad9 elements. The states at the beginning of the time step
are handled by the MGIS material data manager: the Newton loop reverts to
them before each integration and commits them after convergence. The
thermodynamic forces (stresses, including the out-of-plane plane strain
component) and the consistent tangent operator blocks are transferred to the
Gauss point variables of the module, which are then used by the standard
assembly of the module, on the CPU or on an accelerator.

#### Cross validation

The test `2D_Quater_cylinder_von_mises_mgis` runs the quarter cylinder case
with the `VonMisesMGIS` law and compares the displacement field with the
reference values of the native `VonMises` law on the same mesh
([check/2D.quater.cylinder.mgis.t20.txt](check/2D.quater.cylinder.mgis.t20.txt)).
On this case, both implementations agree to machine precision (maximal
relative difference of the order of 1e-13).



## An example elastoplasticity solver
We will provide a simple 2D quarter cylinder example. The physics solved is that the inner boundary of cylinder is subjected to an incremental traction/pressure loading until plastic strains appear. This XML configuration file is used for setting up an elastoplasticity problem simulation in ArcaneFEM. Below is a detailed explanation of each section in the configuration for one of the tests `2D.quater.cylinder.vonMises.arc`.

###### Mesh Configuration

The mesh configuration section specifies the mesh file to be used in the simulation:

```xml
<meshes>
    <mesh>
    <filename>meshes/quater_cylinder.msh</filename>
  </mesh>
</meshes>
```

- **Mesh File:** Defines the mesh file `quater_cylinder.msh` to be used in the simulation. Note that this file should be compatible with version 4.1 `.msh` format from `Gmsh`. The mesh here is a traingular mesh.

###### FEM Configuration

The Finite Element Method (FEM) configuration is provided in the `2D.quater.cylinder.vonMises.arc`.

```xml
<fem>
  <tmax>21.</tmax>
  <dt>1.</dt>
  <constitutive-law>VonMises</constitutive-law>
  <gp-material-tensor-strategy>global</gp-material-tensor-strategy>
  <E>70.0e3</E>
  <nu>0.3</nu>
  <sig0>250.</sig0>
  <f>NULL NULL</f>
  <boundary-conditions>
    <dirichlet>
      <surface>left</surface>
      <value>0.0 NULL</value>
      <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
    </dirichlet>
    <dirichlet>
      <surface>bottom</surface>
      <value>NULL 0.0</value>
      <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
    </dirichlet>
    <traction>
      <surface>inner</surface>
      <traction-input-file>data/traction_quater_cylinder_20steps.txt</traction-input-file>
    </traction>
  </boundary-conditions>
</fem>
```

Let us explain this point wise
- **Incremental loading properties:** The artifitial time step size is `1` and maximun time is `21` for `20` loading steps. 
```xml
  <tmax>21.</tmax>
  <dt>1.</dt>
```
- **Constitutive law:** As of now only `Von Mises` law is implemented and is chosen. 
```xml
  <constitutive-law>VonMises</constitutive-law>
```
- **Quadrature point variable storage strategy:** With `global` we use the ```MRMeshVariables``` of ```Arcane Framework``` to store and reuse state variables on quadrature points. The `local` strategy evaluate these variables in-place at each quadrature point. 
```xml
  <constitutive-law>VonMises</constitutive-law>
  <gp-material-tensor-strategy>global</gp-material-tensor-strategy>
```
- **Material Properties:** The Young's Modulus (E) for the material, defined as `70.0e5`.  The Poisson's Ratio (nu) for the material, defined as `0.3`. The initial state stress, defined as `250`
```xml
<E>21.0e5</E>
<nu>0.28</nu>
<sig0>250.</sig0>
```
- **Source Term / Body Force:** The source term or body force is absent in both x and y directions `NULL NULL`.
```xml
  <f>NULL NULL</f>
```
- **Dirichlet BCs:** The boundary condition on the specified surface (`left`) with given values for `u1` and `u2`, which we set to `0` and `NULL` since the end is clamped in x direction but free to move in y direction. Similarly, The boundary condition on the specified surface (`bottom`) with given values for `u1` and `u2`, which we set to `NULL` and `0` since the end is clamped in y direction but free to move in x direction. Please note that it is expected the mesh `quater_cylinder.msh` contains the physical groups `left` and `bottom`.
```xml
<dirichlet>
  <surface>left</surface>
  <value>0.0 NULL</value>
  <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
</dirichlet>
<dirichlet>
  <surface>bottom</surface>
  <value>NULL 0.0</value>
  <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
</dirichlet>
```
- **Traction BCs:** The incremental load (traction) is applied as an external pressure using data in a file `traction_quater_cylinder_20steps.txt`, in the form of a table. This load is applied on the `inner` boundary of the cylinder representing the curve traced by inner radius. Please note that it is expected the mesh `quater_cylinder.msh` contains the physical group `inner`.
```xml
<traction>
  <surface>inner</surface>
  <traction-input-file>data/traction_quater_cylinder_20steps.txt</traction-input-file>
</traction>
```

###### Post-Processing Configuration

The post-processing configuration is specified to control how and when results are saved:

```xml
<arcane-post-processing>
 <output-period>1</output-period>
 <output>
   <variable>U</variable>
 </output>
</arcane-post-processing>
```

- **Output Period:** Sets the interval at which results are saved.
- **Output Variables:** Lists the variables (`U`) which is the displacement vector to be included in the output.

#### Post Process

For post processing the `ensight.case` file is outputted (in `output/depouillement/` folder), which can be read by ParaView. The output is of the $\mathbb{P}_1$ FE order (on nodes).

#### Cross Validation 
The plot presented here you find a comparative cross-validation of ArcaneFEM elastoplasticity against other reference solvers. 
<img width="500" align="left" src="https://github.com/user-attachments/assets/427c29b1-052a-401c-98b7-1da122eae620" />
