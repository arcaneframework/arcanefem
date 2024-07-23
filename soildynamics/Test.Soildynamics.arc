<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
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
      <filename>bar_dynamic.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>2.</tmax>
    <dt>0.08</dt>
    <alpm>0.20</alpm>
    <alpf>0.40</alpf>
    <rho>1.0</rho>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <penalty>1.e64</penalty>
    <time-discretization>Newmark-beta</time-discretization>
    <result-file>test_soildynamics_results.txt</result-file>
    <dirichlet-boundary-condition>
      <surface>surfaceleft</surface>
      <u1>0.0</u1>
      <u2>0.0</u2>
    </dirichlet-boundary-condition>
    <traction-boundary-condition>
      <surface>surfaceright</surface>
      <t2>0.01</t2>
    </traction-boundary-condition>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <preconditioner>ilu</preconditioner>
    </linear-system>
  </fem>
</case>
