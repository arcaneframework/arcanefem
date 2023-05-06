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
      <filename>ring.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <dirichlet-boundary-condition>
      <surface>inner</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>outer</surface>
      <value>20.0</value>
    </dirichlet-boundary-condition>
  </fem>
</case>
