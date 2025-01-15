<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Cut circle 2D</title>
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
    <result-file>check/poisson_test_ref_circle_2D.txt</result-file>
    <f>5.5</f>
    <boundary-conditions>
      <dirichlet>
        <surface>horizontal</surface>
        <value>0.5</value>
      </dirichlet>
    </boundary-conditions>
    <linear-system name="HypreLinearSystem">
      <rtol>0.</rtol>
      <atol>1e-15</atol>
      <amg-threshold>0.25</amg-threshold>
    </linear-system>
    <bsr>true</bsr>
  </fem>
</case>

