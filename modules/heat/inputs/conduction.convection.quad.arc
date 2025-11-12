<?xml version="1.0"?>
<case codename="Heat" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>HeatLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>2</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/plate.quad.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <lambda>1.75</lambda>
    <tmax>20.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e31</penalty>
        <surface>left</surface>
        <value>10.0</value>
      </dirichlet>
    </boundary-conditions>
    <convection-boundary-condition>
      <surface>right</surface>
      <h>5.</h>
      <Text>15.</Text>
    </convection-boundary-condition>
  </fem>
</case>
