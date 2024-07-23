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
      <filename>plancher.10k.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <dirichlet-point-condition>
      <node>topLeftCorner</node>
      <value>50.0</value>
    </dirichlet-point-condition>
    <dirichlet-point-condition>
      <node>topRightCorner</node>
      <value>20.0</value>
    </dirichlet-point-condition>
    <dirichlet-point-condition>
      <node>botLeftCorner</node>
      <value>20.0</value>
    </dirichlet-point-condition>
    <dirichlet-point-condition>
      <node>botRightCorner</node>
      <value>50.0</value>
    </dirichlet-point-condition>
  </fem>
</case>
