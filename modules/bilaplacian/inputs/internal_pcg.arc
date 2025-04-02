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
    <result-file>check/2d_test.txt</result-file>
    <f>-786.25</f>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>145.5</value>
    </dirichlet-boundary-condition>
    <linear-system name="SequentialBasicLinearSystem">
      <epsilon>1.0e-25</epsilon>
      <solver-method>pcg</solver-method>
    </linear-system>
  </fem>
</case>
