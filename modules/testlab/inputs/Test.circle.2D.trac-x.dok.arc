<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Cut circle 2D with CSR sparse matrix format Gpu compatible and node wise technique. The result of this test is compared with poisson_test_ref_circle_2D.txt</title>
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
      <filename>circle_cut.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>poisson_test_ref_circle_x-trac_2D.txt</result-file>
    <f>5.5</f>
    <dirichlet-boundary-condition>
      <surface>horizontal</surface>
      <value>0.5</value>
    </dirichlet-boundary-condition>
    <neumann-boundary-condition>
      <surface>curved</surface>
      <valueX>2.9e4</valueX>
    </neumann-boundary-condition>
  </fem>
</case>
