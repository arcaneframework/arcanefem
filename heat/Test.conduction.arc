<?xml version="1.0"?>
<case codename="Heat" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>HeatLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>100</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>plate.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>1.75</lambda>
    <tmax>10.</tmax>
    <dt>0.1</dt>
    <h>0.25</h>
    <Text>25.0</Text>
    <Tinit>30.0</Tinit>
    <dirichlet-boundary-condition>
      <surface>left</surface>
      <value>10.0</value>
    </dirichlet-boundary-condition>
  </fem>
</case>
