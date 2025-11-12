<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
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
    <solution-comparison-file>poisson_test_ref_L-shape_2D.txt</solution-comparison-file>
    <f>-5.5</f>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.5</value>
    </dirichlet-boundary-condition>
    <bsr>true</bsr>
    <linear-system name="HypreLinearSystem">
      <solver>bicgstab</solver>
      <rtol>0.</rtol>
      <atol>1e-5</atol>
      <amg-threshold>0.25</amg-threshold>
    </linear-system>
  </fem>
</case>
