<?xml version="1.0"?>
<case codename="Elastodynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ElastodynamicsLoop</timeloop>
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
    <alpm>0.20</alpm>
    <alpf>0.40</alpf>
    <rho>1.0</rho>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <time-discretization>Newmark-beta</time-discretization>
    <result-file>check/2D_elastodynamics_bar_transient_traction.txt</result-file>
    <boundary-conditions>
      <dirichlet>
        <surface>surfaceleft</surface>
        <value>0.0 0.0</value>
      </dirichlet>
    </boundary-conditions>
    <traction-boundary-condition>
      <surface>surfaceright</surface>
      <traction-input-file>data/traction_bar_test_1.txt</traction-input-file>
    </traction-boundary-condition>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <preconditioner>ilu</preconditioner>
    </linear-system>
  </fem>
</case>