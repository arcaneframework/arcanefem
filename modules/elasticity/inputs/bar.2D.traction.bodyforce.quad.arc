<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Traction + Body force + Dirichlet</title>
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
      <filename>meshes/plate.quad.msh</filename>
      <subdivider>
        <nb-subdivision>0</nb-subdivision>
      </subdivider>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>3.33 -6.66</f>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 0.0</value>
      </dirichlet>
      <traction>
        <surface>right</surface>
        <value>1.33 2.13</value>
      </traction>
    </boundary-conditions>
  </fem>
</case>
