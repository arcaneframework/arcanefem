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
    <boundary-conditions>
      <dirichlet-point>
        <node>center</node>
        <value>1.8</value>
      </dirichlet-point>
      <neumann>
        <surface>left</surface>
        <value>3.1</value>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
