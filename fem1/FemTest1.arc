<?xml version="1.0"?>
<case codename="FemTest1" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>FemTest1Loop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>plancher.msh</filename>
    </mesh>
  </meshes>

  <Fem1>
    <lambda>1.75</lambda>
  </Fem1>
</case>
