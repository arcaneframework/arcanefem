<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Testlab: 2D Random Mesh with Neumann Boundary Condition (No Dirichlet, No Linear System Specified)</title>
    <timeloop>TestlabLoop</timeloop>
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
      <filename>random.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <solution-comparison-file>test4_results.txt</solution-comparison-file>
    <neumann-boundary-condition>
      <surface>boundary</surface>
      <valueX>2.0</valueX>
      <valueY>5.0</valueY>
    </neumann-boundary-condition>
  </fem>
</case>
