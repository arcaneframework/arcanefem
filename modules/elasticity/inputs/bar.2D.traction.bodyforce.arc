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
      <subdivider>
        <nb-subdivision>0</nb-subdivision>
      </subdivider>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/elasticity_traction_bodyforce_bar_test_ref.txt</result-file>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>3.33 -6.66</f>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 0.0</value>
      </dirichlet>
    </boundary-conditions>
    <traction-boundary-condition>
      <surface>right</surface>
      <t>1.33 2.13</t>
    </traction-boundary-condition>
  </fem>
</case>
