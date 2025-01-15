<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
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
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <penalty>1.e30</penalty>
    <time-discretization>Newmark-beta</time-discretization>
    <paraxial-boundary-condition>
      <surface>lower</surface>
    </paraxial-boundary-condition>
    <traction-boundary-condition>
      <surface>input</surface>
      <t1>0.01</t1>     
      <t2>0.01</t2>       
    </traction-boundary-condition>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <preconditioner>ilu</preconditioner>
    </linear-system>
  </fem>
</case>
