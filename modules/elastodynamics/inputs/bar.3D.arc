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
      <filename>meshes/bar_dynamic_3D.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <time-discretization>Newmark-beta</time-discretization>
    <tmax>.5</tmax>
    <dt>0.08</dt>
    <rho>1.0</rho>
    <f>131.0e2 113.8e6 567.0e8</f>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <boundary-conditions>
      <dirichlet>
        <surface>surfaceleft</surface>
        <value>0.0 0.0 0.0</value>
      </dirichlet>
      <traction>
        <surface>surfaceright</surface>
        <value>NULL 1869.1e2 NULL</value>
      </traction>
    </boundary-conditions>
  </fem>
</case>
