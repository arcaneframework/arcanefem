<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>2D elasticity with PETSc GAMG and rigid-body near null space</title>
    <timeloop>ElasticityLoop</timeloop>
  </arcane>
  <meshes>
    <mesh>
      <filename>meshes/bar.msh</filename>
    </mesh>
  </meshes>
  <fem>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>NULL -1.0</f>
    <matrix-format>BSR</matrix-format>
    <petsc-flags>-ksp_monitor -mat_type aij</petsc-flags>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 0.0</value>
        <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
      </dirichlet>
    </boundary-conditions>
    <linear-system name="PetscLinearSystem">
      <solver>cg</solver>
      <pc-type>gamg</pc-type>
      <amg-near-null-space>true</amg-near-null-space>
      <rtol>1.0e-9</rtol>
    </linear-system>
  </fem>
</case>
