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
      <filename>meshes/bar_dynamic_3Dhexa.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <time-discretization>Newmark-beta</time-discretization>
    <tmax>.5</tmax>
    <dt>0.08</dt>
    <rho>1232434.0</rho>
    <f>0. 0. 0.</f>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 0.0 0.0</value>
      </dirichlet>
      <traction>
        <surface>right</surface>
        <value>1.7e6, 2.7e6, 2.4e7</value>
      </traction>
    </boundary-conditions>
  </fem>
</case>