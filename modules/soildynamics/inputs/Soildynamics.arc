<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Simple dirichlet traction case</title>
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
      <filename>meshes/bar_dynamic.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>2.</tmax>
    <dt>0.08</dt>
    <rho>1.0</rho>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <result-file>check/test_soildynamics_results.txt</result-file>
    <boundary-conditions>
      <dirichlet>
        <surface>surfaceleft</surface>
        <value>0.0 0.0</value>
      </dirichlet>
    </boundary-conditions>
    <traction-boundary-condition>
      <surface>surfaceright</surface>
      <t>NULL 0.01</t>
    </traction-boundary-condition>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <preconditioner>ilu</preconditioner>
    </linear-system>
  </fem>
</case>
