<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>2D bar bending under own load</title>
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
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>NULL -1.0</f>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 0.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
