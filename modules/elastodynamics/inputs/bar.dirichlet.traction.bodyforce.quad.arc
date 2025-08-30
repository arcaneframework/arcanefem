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
    <result-file>check/2D_elastodynamics_bar_transient_traction.quad.txt</result-file>
    <tmax>1.</tmax>
    <dt>0.08</dt>
    <alpm>0.00</alpm>
    <alpf>0.00</alpf>
    <rho>1.0</rho>
    <f>NULL -2000.</f>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <time-discretization>Newmark-beta</time-discretization>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <surface>surfaceleft</surface>
        <value>0.0 0.0</value>
      </dirichlet>
      <traction>
        <surface>surfaceright</surface>
        <value>NULL -1.</value>
      </traction>
    </boundary-conditions>
    <linear-system>
      <solver-backend>hypre</solver-backend>
    </linear-system>
  </fem>
</case>