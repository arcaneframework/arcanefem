<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Traction + Body force + Dirichlet</title>
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
      <filename>meshes/bar.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/elasticity_traction_bodyforce_bar_test_ref.txt</result-file>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>3.33 -6.66</f>
    <dirichlet-boundary-condition>
      <surface>left</surface>
      <u>0.0 0.0</u>
    </dirichlet-boundary-condition>
    <traction-boundary-condition>
      <surface>right</surface>
      <t1>1.33</t1>
      <t2>2.13</t2>
    </traction-boundary-condition>
  </fem>
</case>
