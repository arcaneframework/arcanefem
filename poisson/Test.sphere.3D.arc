<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PoissonLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>sphere_cut.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>test_sphere_3D_results.txt</result-file>
    <blcsr>true</blcsr>
    <f>-0.01</f>
    <mesh-type>TETRA4</mesh-type>
    <dirichlet-boundary-condition>
      <surface>sphere</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Cut</surface>
      <value>10.0</value>
    </dirichlet-boundary-condition>
    <linear-system name="HypreLinearSystem"/>
  </fem>
</case>
