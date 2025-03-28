<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3D soildynamics test</title>
    <timeloop>SoildynamicsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
     <variable>V</variable>
     <variable>A</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/cube_double_couple_3d.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>.1</tmax>
    <dt>0.01</dt>
    <cs>21.5</cs>
    <cp>40.</cp>
    <rho>9.0</rho>
    <paraxial-boundary-condition>
      <surface>leftsur</surface>
    </paraxial-boundary-condition>
    <double-couple>
      <north-node-name>dcNorth</north-node-name>
      <south-node-name>dcSouth</south-node-name>
      <east-node-name>dcEast</east-node-name>
      <west-node-name>dcWest</west-node-name>
      <method>force-based</method>
      <double-couple-input-file>data/force_loading_dc.txt</double-couple-input-file>
    </double-couple>
    <result-file>check/3d_test_paraxial_double_couple.txt</result-file>
  </fem>
</case>
