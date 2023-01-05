<?xml version='1.0'?>
<case codename="Passmo" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PassmoLoop</timeloop>
  </arcane>
  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>Displacement</variable>
   </output>
  </arcane-post-processing>

    <mesh>
      <filename>plancher.msh</filename>
    </mesh>

    <elastodynamic>
      <boundary-condition>
        <surface>XMAX</surface>
        <value>50.0</value>
      </boundary-condition>
    </elastodynamic>
</case>
