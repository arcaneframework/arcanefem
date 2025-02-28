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
      <filename>meshes/sphere_cut.msh</filename>
      <subdivider>
        <nb-subdivision>0</nb-subdivision>
      </subdivider>
    </mesh>
  </meshes>

  <fem>
    <E>21.0e3</E>
    <nu>0.28</nu>
    <f>-19 -23 -42</f>
    <dirichlet-boundary-condition>
      <surface>verticalXY</surface>
      <u>NULL NULL 0.0</u>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>verticalYZ</surface>
      <u>0.0 NULL NULL</u>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>horizontal</surface>
      <u>NULL 0.0 NULL</u>
    </dirichlet-boundary-condition>
    <traction-boundary-condition>
      <surface>curved</surface>
      <t>1.99 1.5 2.4</t>
    </traction-boundary-condition>
  </fem>
</case>
