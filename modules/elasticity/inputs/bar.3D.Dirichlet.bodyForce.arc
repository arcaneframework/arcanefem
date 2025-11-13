<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3D Linear Elasticity</title>
    <timeloop>ElasticityLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/bar_dynamic_3D.msh</filename>
      <subdivider>
        <nb-subdivision>0</nb-subdivision>
      </subdivider>
    </mesh>
  </meshes>

  <fem>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>-1.0</f>
    <boundary-conditions>
      <dirichlet>
        <surface>surfaceleft</surface>
        <value>0.0 0.0 0.0</value>
      </dirichlet>
      <dirichlet>
        <surface>surfaceright</surface>
        <value>NULL 1.0 NULL</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
