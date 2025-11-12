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
      <filename>meshes/semi-circle-soil.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>0.08</tmax>
    <dt>0.01</dt>
    <E>6.62e6</E>
    <nu>0.45</nu>
    <rho>2500.0</rho>
    <f>3359.6 3452.3</f>
    <paraxial-boundary-condition>
      <surface>lower</surface>
    </paraxial-boundary-condition>
    <boundary-conditions>
      <dirichlet-point>
        <node>source</node>
        <value>0.0 0.0003</value>
      </dirichlet-point>
      <traction>
        <surface>input</surface>
        <value>0.01 0.01</value>
      </traction>
    </boundary-conditions>
  </fem>
</case>
