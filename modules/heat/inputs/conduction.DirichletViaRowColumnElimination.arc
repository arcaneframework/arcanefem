<?xml version="1.0"?>
<case codename="Heat" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>HeatLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>2</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/plate.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/2d_conduction.txt</result-file>
    <lambda>1.75</lambda>
    <tmax>20.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
        <surface>left</surface>
        <value>10.0</value>
      </dirichlet>
    </boundary-conditions>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <solver-method>pcg</solver-method>
      <preconditioner>amg</preconditioner>
    </linear-system>
  </fem>
</case>
