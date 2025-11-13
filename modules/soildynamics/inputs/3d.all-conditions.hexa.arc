<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3D dummy test to pass via all bc</title>
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
      <filename>meshes/cube_double_couple_3d.hexa.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>.1</tmax>
    <dt>0.01</dt>
    <cs>21.5</cs>
    <cp>40.</cp>
    <rho>9.0</rho>
    <f>0.03 0.04 18.0</f>
    <hex-quad-mesh>true</hex-quad-mesh>
    <boundary-conditions>
      <dirichlet>
        <surface>boundaries</surface>
        <value>0.0 0.0 0.0</value>
      </dirichlet>
      <dirichlet-point>
        <node>north</node>
        <value>0.0 NULL NULL</value>
      </dirichlet-point>
      <traction>
        <surface>boundaries</surface>
        <value>0.01 0.01 0.02</value>
      </traction>
    </boundary-conditions>
    <paraxial-boundary-condition>
      <surface>boundaries</surface>
    </paraxial-boundary-condition>
    <double-couple>
      <north-node-name>north</north-node-name>
      <south-node-name>south</south-node-name>
      <east-node-name>east</east-node-name>
      <west-node-name>west</west-node-name>
      <method>force-based</method>
      <double-couple-input-file>data/force_loading_dc.txt</double-couple-input-file>
    </double-couple>
    <linear-system>
      <solver-backend>hypre</solver-backend>
    </linear-system>
  </fem>
</case>
