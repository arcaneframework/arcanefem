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
     <variable>u1</variable>
     <variable>u2</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/bilap.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>-1.0</f>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.05</value>
    </dirichlet-boundary-condition>
    <linear-system name="SequentialBasicLinearSystem">
      <epsilon>1.0e-25</epsilon>
      <solver-method>pcg</solver-method>
    </linear-system>
  </fem>
</case>
