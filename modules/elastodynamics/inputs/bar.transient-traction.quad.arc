<?xml version="1.0"?>
<case codename="Elastodynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ElastodynamicsLoop</timeloop>
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
      <filename>meshes/bar_dynamic_quad.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <result-file>check/bar_transient-traction.quad.txt</result-file>
    <tmax>2.</tmax>
    <dt>0.08</dt>
    <alpm>0.20</alpm>
    <alpf>0.40</alpf>
    <rho>1.0</rho>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <time-discretization>Newmark-beta</time-discretization>
    <boundary-conditions>
      <dirichlet>
        <surface>surfaceleft</surface>
        <value>0.0 0.0</value>
      </dirichlet>
      <traction>
        <surface>surfaceright</surface>
        <traction-input-file>data/traction_bar_test_1.txt</traction-input-file>
      </traction>
    </boundary-conditions>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <preconditioner>ic</preconditioner>
    </linear-system>
  </fem>
</case>