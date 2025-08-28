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
    <time-discretization>Newmark-beta</time-discretization>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <surface>bottom</surface>
        <value>0.0 0.0 0.0</value>
      </dirichlet>
      <dirichlet-point>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <node>center</node>
        <value>18.0 13.0 14.0</value>
      </dirichlet-point>
    </boundary-conditions>
  </fem>
</case>