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
    <f>0.0 0.0 0.0</f>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <penalty>1.e64</penalty>
    <dirichlet-boundary-condition>
      <surface>surfaceleft</surface>
      <u>0.0 0.0 0.0</u>
    </dirichlet-boundary-condition>
    <traction-boundary-condition>
      <surface>surfaceright</surface>
      <t>NULL 0.01 NULL</t>
    </traction-boundary-condition>
  </fem>
</case>
