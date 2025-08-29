<?xml version="1.0"?>
<case codename="Elastodynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>DIrichlet and source term for quad mesh</title>
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
    <time-discretization>Newmark-beta</time-discretization>
    <result-file>check/bar.quad.txt</result-file>
    <tmax>2.</tmax>
    <dt>0.08</dt>
    <rho>12.0</rho>
    <f>0. 13.5e2</f>
    <lambda>576.9230769</lambda>
    <mu>384.6153846</mu>
    <boundary-conditions>
      <dirichlet>
        <surface>surfaceleft</surface>
        <value>0.0 0.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>