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
      <filename>meshes/truncated_cube.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/test_3D_truncated-cube.txt</result-file>
    <boundary-conditions>
      <dirichlet-point>
        <node>center</node>
        <value>18.8</value>
      </dirichlet-point>
      <dirichlet>
        <surface>left</surface>
        <value>0.8</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
