<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PoissonLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>random.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>test4_results.txt</result-file>
    <neumann-boundary-condition>
      <surface>boundary</surface>
      <valueX>2.0</valueX>
      <valueY>5.0</valueY>
    </neumann-boundary-condition>
  </fem>
</case>
