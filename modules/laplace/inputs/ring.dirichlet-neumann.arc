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
      <filename>meshes/ring.msh</filename>
    </mesh>
  </meshes>

  <fem>
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
