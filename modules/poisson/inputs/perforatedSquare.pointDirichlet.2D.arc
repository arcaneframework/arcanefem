<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Poisson Solver With Point Dirichlet Conditions</title>
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
      <filename>meshes/plancher.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/poisson_test_point_dirichlet_2D.txt</result-file>
    <boundary-conditions>
      <dirichlet-point>
        <node>topLeftCorner</node>
        <value>50.0</value>
      </dirichlet-point>
      <dirichlet-point>
        <node>topRightCorner</node>
        <value>20.0</value>
      </dirichlet-point>
      <dirichlet-point>
        <node>botLeftCorner</node>
        <value>20.0</value>
      </dirichlet-point>
      <dirichlet-point>
        <node>botRightCorner</node>
        <value>50.0</value>
      </dirichlet-point>
    </boundary-conditions>
    <linear-system name="SequentialBasicLinearSystem" />
  </fem>
</case>
