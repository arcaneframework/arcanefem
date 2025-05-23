<?xml version="1.0"?>
<case codename="Electrostatics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ElectrostaticsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>Phi</variable>
     <variable>E</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/box-rod-circle.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <rho>0.0</rho>
    <epsilon>1.0</epsilon>
    <result-file>check/test_2.txt</result-file>

    <boundary-conditions>
      <dirichlet>
        <surface>rod1</surface>
        <value>-1.0</value>
      </dirichlet>
      <dirichlet>
        <surface>circle</surface>
        <value>1.0</value>
      </dirichlet>
      <dirichlet>
        <surface>external</surface>
        <value>0.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
