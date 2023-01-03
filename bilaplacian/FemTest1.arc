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
     <variable>NodeTemp</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>bilap.msh</filename>
    </mesh>
  </meshes>

  <fem1>
    <lambda>1.</lambda>
    <qdot>-1.0</qdot>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.05</value>
    </dirichlet-boundary-condition>
  </fem1>
</case>
