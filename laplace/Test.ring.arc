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
      <filename>ring.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <boundary-conditions>
      <dirichlet>
        <surface>inner</surface>
        <value>50.0</value>
      </dirichlet>
      <dirichlet>
        <surface>outer</surface>
        <value>20.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
