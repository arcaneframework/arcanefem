<?xml version="1.0"?>
<case codename="Heat" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>HeatLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>2</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/plate.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/2d_conduction_convection.txt</result-file>
    <lambda>1.75</lambda>
    <tmax>20.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
    <dirichlet-boundary-condition>
      <surface>left</surface>
      <value>10.0</value>
    </dirichlet-boundary-condition>
    <convection-boundary-condition>
      <surface>right</surface>
      <h>1.</h>
      <Text>20.</Text>
    </convection-boundary-condition>
    <convection-boundary-condition>
      <surface>top</surface>
      <h>1.</h>
      <Text>20.</Text>
    </convection-boundary-condition>
    <convection-boundary-condition>
      <surface>bottom</surface>
      <h>1.</h>
      <Text>20.</Text>
    </convection-boundary-condition>
  </fem>
</case>
