<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Double couple test with paraxial boundries</title>
    <timeloop>SoildynamicsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
     <variable>V</variable>
     <variable>A</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/soil_2d.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>0.2</tmax>
    <dt>0.01</dt>
    <cs>11.5</cs>
    <cp>20.</cp>
    <rho>9.0</rho>
    <double-couple>
      <north-node-name>DcNorthPointCord</north-node-name>
      <south-node-name>DcSouthPointCord</south-node-name>
      <east-node-name>DcEastPointCord</east-node-name>
      <west-node-name>DcWestPointCord</west-node-name>
      <method>force-based</method>
      <double-couple-input-file>data/force_loading_dc.txt</double-couple-input-file>
    </double-couple>
    <paraxial-boundary-condition>
      <surface>surfaceLeft</surface>
    </paraxial-boundary-condition>
    <paraxial-boundary-condition>
      <surface>surfaceRight</surface>
    </paraxial-boundary-condition>
    <paraxial-boundary-condition>
      <surface>surfaceBottom</surface>
    </paraxial-boundary-condition>
    <linear-system>
      <solver-backend>hypre</solver-backend>
    </linear-system>
  </fem>
</case>
