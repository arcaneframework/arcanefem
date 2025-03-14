<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3D sdummy test to pass via all bc</title>
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
    <f>0.03 0.04 18.0</f>
    <traction-boundary-condition>
      <surface>leftsur</surface>
      <t>0.01 0.01 0.02</t>
    </traction-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>rightsur</surface>
      <u>0.0 0.0 0.0</u>
    </dirichlet-boundary-condition>
    <traction-boundary-condition>
      <surface>topsur</surface>
      <traction-input-file>data/3d-traction.txt</traction-input-file>
    </traction-boundary-condition>
    <paraxial-boundary-condition>
      <surface>botsur</surface>
    </paraxial-boundary-condition>
    <dirichlet-point-condition>
      <node>dctop</node>
      <u>0.0 NULL NULL</u>
    </dirichlet-point-condition>
    <double-couple>
      <north-node-name>dcNorth</north-node-name>
      <south-node-name>dcSouth</south-node-name>
      <east-node-name>dcEast</east-node-name>
      <west-node-name>dcWest</west-node-name>
      <method>force-based</method>
      <double-couple-input-file>data/force_loading_dc.txt</double-couple-input-file>
    </double-couple>
    <linear-system>
      <solver-backend>hypre</solver-backend>
    </linear-system>
  </fem>
</case>
