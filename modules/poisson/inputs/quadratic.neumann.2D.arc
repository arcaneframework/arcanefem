<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Quadratic quadrilateral Neumann validation</title>
    <timeloop>PoissonLoop</timeloop>
  </arcane>

  <meshes>
    <mesh>
      <filename>meshes/mms.quad8_n4.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <!-- Exact solution u=x^2: -Delta(u)=-2, u=0 at x=0,
         and the outward scalar flux is 2 at x=1. -->
    <matrix-format>DOK</matrix-format>
    <f>-2.0</f>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0</value>
      </dirichlet>
      <neumann>
        <surface>right</surface>
        <value>2.0</value>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
