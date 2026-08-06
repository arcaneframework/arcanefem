<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Cube 3D hexahedral elements</title>
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
      <filename>meshes/3x3x3_cube_hexa8.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <f>9.8</f>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.5</value>
      </dirichlet>
      <neumann>
        <surface>right</surface>
        <value>13.9</value>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
