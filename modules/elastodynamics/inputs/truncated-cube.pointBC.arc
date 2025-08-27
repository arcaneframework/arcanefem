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
      <filename>meshes/truncated_cube.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/truncated-cube_point-bc.txt</result-file>
    <tmax>1.</tmax>
    <dt>0.08</dt>
    <rho>1.0</rho>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <f>145.5e5 56456.5e6 87842.5e5</f>
    <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
    <time-discretization>Newmark-beta</time-discretization>
    <dirichlet-boundary-condition>
      <surface>bottom</surface>
      <u>0.0 0.0 0.0</u>
    </dirichlet-boundary-condition>
    <dirichlet-point-condition>
      <node>center</node>
      <u>18.0 13.0 14.0</u>
    </dirichlet-point-condition>
    <linear-system>
      <solver-backend>hypre</solver-backend>
      <solver-method>bicgstab</solver-method>
      <epsilon>1e-8</epsilon>
    </linear-system>
  </fem>
</case>