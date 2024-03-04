<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>SoildynamicsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>5</output-period>
   <output>
     <variable>U</variable>
     <variable>V</variable>
     <variable>A</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>square_double-couple.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>2.0</tmax>
    <dt>0.01</dt>
    <cs>2</cs>
    <cp>4</cp>
    <rho>1</rho>
    <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
    <time-discretization>Newmark-beta</time-discretization>
    <result-file>test_paraxial_results.txt</result-file>
    <double-couple>
      <north-node-name>sourceT</north-node-name>
      <south-node-name>sourceB</south-node-name>
      <east-node-name>sourceR</east-node-name>
      <west-node-name>sourceL</west-node-name>
      <method>force-based</method>
      <double-couple-input-file>force_loading_dc</double-couple-input-file>
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
    <linear-system>
      <solver-backend>hypre</solver-backend>
    </linear-system>
  </fem>
</case>
