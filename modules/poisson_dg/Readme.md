# Solving Poisson equation with Discontinuous Galerkin method ArcaneFEM

<img align="left" width="400" src="https://github.com/arcaneframework/arcanefem/assets/52162083/3646fbb9-5c82-4807-b025-7f3b4c899ca7" alt="poisson_1_large" />

Here, we utilize ArcaneFEM to solve the Poisson equation, which is a fundamental elliptic partial differential equation (PDE). The provided code demonstrates a straightforward implementation of a 2D/3D unstructured mesh Discontinous Galerkin finite element method (DG-FEM) solver on arbitary meshes. Although we shall explain here only 2D for keeping the text simple. 

The Poisson equation is encountered in various physical scenarios, including heat conduction, substance diffusion, membrane elasticity, inviscid fluid flow, electrostatics, twisting of elastic rods, and water waves. It serves as a vital tool for modeling and understanding these phenomena.

## Theory

#### Problem description

The 2D Poisson's equation is solved for a closed meshed domain $\Omega^h$ in order to know the solution $u(x,y)$ within the domain. The equation reads

$$
\frac{\partial}{\partial x}\left( \frac{\partial u}{\partial x} \right) + \frac{\partial}{\partial y}\left( \frac{\partial u}{\partial y} \right) = {\mathcal{f}}   \quad \forall (x,y)\in\Omega^h
$$

or in a more compact form

$$
\nabla^2 u = {\mathcal{f}} \quad \forall (x,y)\in\Omega^h.
$$

To complete the problem description,   first type (Dirichlet) boundary conditions is applied to this problem:

$u = 0.5 \quad \forall(x,y)\in\partial\Omega^h_{\text{horizontal}}\subset\partial \Omega^h,$

Finally the right hand side source is present within the domain

${\mathcal{f}}=5.5$

#### Post Process

<img align="left" width="200" src="https://github.com/user-attachments/assets/66b9449e-e2f7-4607-b910-231def7d2f67" alt="poisson_1_large" />
<img align="left" width="200" src="https://github.com/user-attachments/assets/a68dd3d8-3f9f-424e-8d6a-33e34e41c04b" alt="poisson_1_large" />
Post processing is controlled via

```xml
  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>
```

For post processing the `Mesh0.hdf` file is outputted (in `output/depouillement/vtkhdfv2` folder), which can be read by PARAVIS. The output is of the $\mathbb{P}_1$ FE order (on nodes).
