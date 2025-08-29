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
    <result-file>check/bar_3d.hexa.txt</result-file>
    <time-discretization>Newmark-beta</time-discretization>
    <tmax>.5</tmax>
    <dt>0.08</dt>
    <rho>1232434.0</rho>
    <f>131.0e2 113.8e6 567.0e8</f>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 0.0 0.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
