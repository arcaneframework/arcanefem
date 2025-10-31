<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Constant traction test with point bounday consitions</title>
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
      <filename>meshes/semi-circle-soil-quad.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <tmax>0.2</tmax>
    <dt>0.01</dt>
    <E>6.62e7</E>
    <nu>0.3</nu>
    <rho>250.0</rho>
    <paraxial-boundary-condition>
      <surface>lower</surface>
    </paraxial-boundary-condition>
    <boundary-conditions>
      <dirichlet-point>
        <node>source</node>
        <value>0.0 0.0003</value>
      </dirichlet-point>
    </boundary-conditions>
  </fem>
</case>
