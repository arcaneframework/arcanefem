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
      <filename>meshes/semi-circle.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/semi-ciricle_point-bc.txt</result-file>
    <tmax>1.</tmax>
    <dt>0.08</dt>
    <rho>1.0</rho>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
    <time-discretization>Newmark-beta</time-discretization>
    <dirichlet-boundary-condition>
      <surface>boderCircle</surface>
      <u>0.0 0.0</u>
    </dirichlet-boundary-condition>
    <dirichlet-point-condition>
      <node>source</node>
      <u>10.0 10.0</u>
    </dirichlet-point-condition>
  </fem>
</case>