<?xml version="1.0"?>
<case codename="Elastoplasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Elastoplasticity with  Dirichlet, traction, and bodyforce</title>
    <timeloop>ElastoplasticityLoop</timeloop>
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
    <hex-quad-mesh>true</hex-quad-mesh>
    <E>200e9</E>
    <nu>0.3</nu>
    <f>9.8e9 7.3e9</f>
    <boundary-conditions>
      <dirichlet>
        <surface>bot</surface>
        <value>0.0 0.0</value>
      </dirichlet>
      <traction>
        <surface>top</surface>
        <value>13.3e9 14.5e5</value>
      </traction>
    </boundary-conditions>
  </fem>
</case>