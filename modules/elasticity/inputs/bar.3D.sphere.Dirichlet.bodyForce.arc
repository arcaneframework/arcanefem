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
    <boundary-conditions>
      <dirichlet>
        <surface>verticalXY</surface>
        <value>NULL NULL 0.0</value>
      </dirichlet>
      <dirichlet>
        <surface>verticalYZ</surface>
        <value>0.0 NULL NULL</value>
      </dirichlet>
      <dirichlet>
        <surface>horizontal</surface>
        <value>NULL 0.0 NULL</value>
      </dirichlet>
    </boundary-conditions>
    <traction-boundary-condition>
      <surface>curved</surface>
      <t>1.99 1.5 2.4</t>
    </traction-boundary-condition>
  </fem>
</case>
