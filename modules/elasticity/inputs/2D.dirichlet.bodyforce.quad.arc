<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Elasticity with two Dirichlet and bodyforce</title>
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
      <filename>meshes/five_quads.msh</filename>
      <subdivider>
        <nb-subdivision>0</nb-subdivision>
      </subdivider>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/dirichlet_bodyforce.quad.txt</result-file>
    <hex-quad-mesh>true</hex-quad-mesh>
    <E>200e9</E>
    <nu>0.3</nu>
    <f>-9818949214245.0, -7818949234281.0</f>
    <boundary-conditions>
      <dirichlet>
        <surface>bot</surface>
        <value>0.0 0.0</value>
      </dirichlet>
      <dirichlet>
        <surface>top</surface>
        <value>1.9 14.5</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>