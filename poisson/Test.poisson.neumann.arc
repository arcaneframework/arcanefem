<?xml version="1.0"?>
<case codename="FemTest" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>FemTestLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>random.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>1.0</lambda>
    <neumann-boundary-condition>
      <surface>boundary</surface>
      <valueX>2.0</valueX>
      <valueY>5.0</valueY>
    </neumann-boundary-condition>
  </fem>
</case>
