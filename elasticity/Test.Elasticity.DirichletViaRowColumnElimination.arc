<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
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
      <filename>bar.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f2>-1.0</f2>
    <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
    <dirichlet-boundary-condition>
      <surface>left</surface>
      <u1>0.0</u1>
      <u2>0.0</u2>
    </dirichlet-boundary-condition>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <solver-method>cg</solver-method>
      <preconditioner>amg</preconditioner>
    </linear-system>
  </fem>
</case>
