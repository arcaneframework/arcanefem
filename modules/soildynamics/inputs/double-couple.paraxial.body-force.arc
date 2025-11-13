<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Double couple test with paraxial boundries and body force</title>
    <timeloop>SoildynamicsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>5</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
     <variable>V</variable>
     <variable>A</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/square_double-couple.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>0.2</tmax>
    <dt>0.01</dt>
    <cs>2</cs>
    <cp>4</cp>
    <rho>1</rho>
    <f>1255.1 32289.5</f>
    <double-couple>
      <north-node-name>sourceT</north-node-name>
      <south-node-name>sourceB</south-node-name>
      <east-node-name>sourceR</east-node-name>
      <west-node-name>sourceL</west-node-name>
      <method>force-based</method>
      <double-couple-input-file>data/force_loading_dc.txt</double-couple-input-file>
    </double-couple>
    <paraxial-boundary-condition>
      <surface>left</surface>
    </paraxial-boundary-condition>
    <paraxial-boundary-condition>
      <surface>top</surface>
    </paraxial-boundary-condition>
    <paraxial-boundary-condition>
      <surface>right</surface>
    </paraxial-boundary-condition>
    <paraxial-boundary-condition>
      <surface>bottom</surface>
    </paraxial-boundary-condition>
  </fem>
</case>
