<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>L-Shape 2D with COO sparse matrix format GPU compatible. The result of this test is compared with poisson_test_ref_L-shape_2D.txt</title>
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
      <filename>L-shape.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <coo-gpu>true</coo-gpu>
    <solution-comparison-file>poisson_test_ref_L-shape_2D.txt</solution-comparison-file>
    <f>-5.5</f>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.5</value>
    </dirichlet-boundary-condition>
  </fem>
</case>
