<?xml version="1.0"?>
<case codename="Bilaplacian" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>BilaplacianLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>u1</variable>
     <variable>u2</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>bilap.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>-1.0</f>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.05</value>
    </dirichlet-boundary-condition>
    <linear-system name="SequentialBasicLinearSystem">
      <solver-method>direct</solver-method>
    </linear-system>
  </fem>
</case>
