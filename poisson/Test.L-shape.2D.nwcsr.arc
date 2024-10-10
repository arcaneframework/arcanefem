<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>L-shape 2D with CSR sparse matrix format Gpu compatible and node wise technique. The result of this test is compared with test_poisson_results.txt</title>
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
      <filename>L-shape.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <nwcsr>true</nwcsr>
    <result-file>test_poisson_results.txt</result-file>
    <f>-1.0</f>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
  </fem>
</case>
