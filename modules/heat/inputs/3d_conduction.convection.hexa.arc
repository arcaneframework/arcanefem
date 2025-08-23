<?xml version="1.0"?>
<case codename="Heat" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3D heat conduction poroblem</title>
    <timeloop>HeatLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/truncated_cube.hexa.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/3d_conduction_convection.hexa.txt</result-file>
    <hex-quad-mesh>true</hex-quad-mesh>
    <lambda>1.75</lambda>
    <tmax>20.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>10.</value>
      </dirichlet>
    </boundary-conditions>
    <convection-boundary-condition>
      <surface>right</surface>
      <h>5.</h>
      <Text>15.</Text>
    </convection-boundary-condition>
  </fem>
</case>