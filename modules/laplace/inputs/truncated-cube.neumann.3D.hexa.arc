<?xml version="1.0"?>
<case codename="Laplace" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>LaplaceLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/truncated_cube.hexa.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <boundary-conditions>
      <dirichlet>
        <surface>horizontal</surface>
        <value>1.8</value>
      </dirichlet>
      <neumann>
        <surface>bottom</surface>
        <value>3.1</value>
      </neumann>
    </boundary-conditions>
    <result-file>check/test_trucated-cube_hexa.txt</result-file>
  </fem>
</case>
