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
    </mesh>
  </meshes>

  <fem>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>-1.0</f>
    <dirichlet-boundary-condition>
      <surface>surfaceleft</surface>
      <u>0.0 0.0 0.0</u>
    </dirichlet-boundary-condition>
    <traction-boundary-condition>
      <surface>surfaceright</surface>
      <t>NULL 1.0 NULL</t>
    </traction-boundary-condition>
    <enforce-Dirichlet-method>RowElimination</enforce-Dirichlet-method>
    <linear-system name="HypreLinearSystem">
      <solver>bicgstab</solver>
      <rtol>1e-9</rtol>
      <atol>0.</atol>
    </linear-system>
  </fem>
</case>
