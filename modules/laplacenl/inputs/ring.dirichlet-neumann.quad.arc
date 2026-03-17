<?xml version="1.0"?>
<case codename="Laplace" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Ring geomerty with quad mesh and Dirichlet and Neumann conditions</title>
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
      <filename>meshes/ring.quad.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <boundary-conditions>
      <dirichlet>
        <surface>inner</surface>
        <value>50.0</value>
      </dirichlet>
      <neumann>
        <surface>outer</surface>
        <value>17.8</value>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
