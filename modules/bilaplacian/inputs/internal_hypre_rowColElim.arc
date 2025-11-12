<?xml version="1.0"?>
<case codename="Bilaplacian" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>BilaplacianLoop</timeloop>
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
      <filename>meshes/bilap.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>-786.25</f>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <u>145.5</u>
    </dirichlet-boundary-condition>
      <neumann-boundary-condition>
      <surface>boundary</surface>
      <value>0</value>
    </neumann-boundary-condition>
    <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
    <linear-system>
      <solver-backend>hypre</solver-backend>
      <solver-method>gmres</solver-method>
      <epsilon>1e-15</epsilon>
    </linear-system>
  </fem>
</case>
