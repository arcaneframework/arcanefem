<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
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
    <result-file>check/elasticity_point-dirichlet_bar_test_ref.txt</result-file>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>NULL -1.0</f>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 NULL</value>
      </dirichlet>
      <dirichlet>
        <surface>right</surface>
        <value>1.0 NULL</value>
      </dirichlet>
      <dirichlet-point>
        <node>botLeft</node>
        <value>0.0 0.0</value>
      </dirichlet-point>
      <dirichlet-point>
        <node>botRight</node>
        <value>NULL 0.0</value>
      </dirichlet-point>
    </boundary-conditions>
  </fem>
</case>
