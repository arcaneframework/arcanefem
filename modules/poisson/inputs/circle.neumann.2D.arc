<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Cut circle 2D with Dirichlet and Neumann</title>
    <timeloop>PoissonLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/circle_cut.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/poisson_test_ref_circle_neumann_2D.txt</result-file>
    <f>5.5</f>
    <boundary-conditions>
      <dirichlet>
        <surface>horizontal</surface>
        <value>0.5</value>
      </dirichlet>
      <neumann>
        <surface>curved</surface>
        <valueX>-0.35</valueX>
        <valueY>1.65</valueY>
      </neumann>
    </boundary-conditions>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <epsilon>1e-15</epsilon>
    </linear-system>
  </fem>
</case>
