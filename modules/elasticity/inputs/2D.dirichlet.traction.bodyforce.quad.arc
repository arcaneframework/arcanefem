<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Elasticity with  Dirichlet, traction, and bodyforce</title>
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
    <result-file>check/dirichlet_traction_bodyforce.quad.txt</result-file>
    <hex-quad-mesh>true</hex-quad-mesh>
    <E>200e9</E>
    <nu>0.3</nu>
    <f>9.8e9 7.3e9</f>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <dirichlet-boundary-condition>
      <surface>bot</surface>
      <u>0.0 0.0</u>
    </dirichlet-boundary-condition>
    <traction-boundary-condition>
      <surface>top</surface>
      <t>13.3e9 14.5e5</t>
    </traction-boundary-condition>
  </fem>
</case>