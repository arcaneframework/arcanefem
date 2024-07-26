<?xml version="1.0"?>
<case codename="Acoustics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>AcousticsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>sub.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <kc2>.11e1</kc2>
    <neumann-boundary-condition>
      <surface>inner1</surface>
      <value>1.0</value>
    </neumann-boundary-condition>
    <linear-system name="SequentialBasicLinearSystem" />
  </fem>
</case>
