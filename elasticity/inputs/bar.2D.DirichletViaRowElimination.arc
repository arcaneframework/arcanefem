<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ElasticityLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/bar.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>NULL -1.0</f>
    <enforce-Dirichlet-method>RowElimination</enforce-Dirichlet-method>
    <dirichlet-boundary-condition>
      <surface>left</surface>
      <u>0.0 0.0</u>
    </dirichlet-boundary-condition>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <solver-method>gmres</solver-method>
      <preconditioner>ilu</preconditioner>
    </linear-system>
  </fem>
</case>