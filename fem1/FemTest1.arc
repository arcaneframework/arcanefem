<?xml version="1.0"?>
<case codename="FemTest1" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>FemTest1Loop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>plancher.msh</filename>
    </mesh>
  </meshes>

  <fem1>
    <lambda>1.75</lambda>
    <dirichlet-boundary-condition>
      <surface>Cercle</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Bas</surface>
      <value>5.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Haut</surface>
      <value>21.0</value>
    </dirichlet-boundary-condition>
  </fem1>
</case>
