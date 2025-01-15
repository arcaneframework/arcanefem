<?xml version="1.0"?>
<case codename="Laplace" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>LaplaceLoop</timeloop>
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
    <result-file>check/test3_results.txt</result-file>
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
    <bsr>true</bsr>
    <linear-system name="HypreLinearSystem">
      <rtol>0.</rtol>
      <atol>1e-15</atol>
      <amg-threshold>0.25</amg-threshold>
    </linear-system>
  </fem>
</case>
