<?xml version="1.0"?>
<case codename="Electrostatics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ElectrostaticsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>L-shape.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <rho>-1.0</rho>
    <epsilon>1.0</epsilon>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
  </fem>
</case>
