<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3D elasticity with Hypre BoomerAMG and rigid-body near null space</title>
    <timeloop>ElasticityLoop</timeloop>
  </arcane>
  <meshes>
    <mesh>
      <filename>meshes/bar_dynamic_3D.msh</filename>
    </mesh>
  </meshes>
  <fem>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>-1.0</f>
    <matrix-format>BSR</matrix-format>
    <boundary-conditions>
      <dirichlet>
        <surface>surfaceleft</surface>
        <value>0.0 0.0 0.0</value>
        <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
      </dirichlet>
      <dirichlet>
        <surface>surfaceright</surface>
        <value>NULL 1.0 NULL</value>
        <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
      </dirichlet>
    </boundary-conditions>
    <linear-system name="HypreLinearSystem">
      <solver>cg</solver>
      <preconditioner>amg</preconditioner>
      <amg-near-null-space>true</amg-near-null-space>
      <rtol>1.0e-9</rtol>
    </linear-system>
  </fem>
</case>
