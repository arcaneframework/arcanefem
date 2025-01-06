<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Traction + Body force + Dirichlet using BSR Format</title>
    <timeloop>ElasticityLoop</timeloop>
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
      <filename>bar.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>elasticity_traction_bodyforce_bar_test_ref.txt</result-file>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f1>3.33</f1>
    <f2>-6.66</f2>
    <dirichlet-boundary-condition>
      <surface>left</surface>
      <u1>0.1</u1>
      <u2>0.1</u2>
    </dirichlet-boundary-condition>
    <traction-boundary-condition>
      <surface>right</surface>
      <t1>1.33</t1>
    </traction-boundary-condition>
    <bsr>true</bsr>
  </fem>
</case>
