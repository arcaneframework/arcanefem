<?xml version="1.0"?>
<case codename="Heat" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3D heat conduction poroblem</title>
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
      <filename>meshes/truncated_cube.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>1.75</lambda>
    <tmax>1.</tmax>
    <dt>0.1</dt>
    <Tinit>30.0</Tinit>
    <result-file>check/3d_conduction_convection_pointBC.txt</result-file>
    <convection-boundary-condition>
      <surface>top</surface>
      <h>1.9</h>
      <Text>30.</Text>
    </convection-boundary-condition>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <surface>bottom</surface>
        <value>2.6</value>
      </dirichlet>
      <dirichlet-point>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <node>center</node>
        <value>100.6</value>
      </dirichlet-point>
    </boundary-conditions>
  </fem>
</case>
