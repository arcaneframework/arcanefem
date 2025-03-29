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
      <filename>meshes/multi-material.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>20.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
    <material-property>
      <volume>Mat1</volume>
      <lambda>1.3</lambda>
    </material-property> 
    <material-property>
      <volume>Mat2</volume>
      <lambda>12.6</lambda>
    </material-property> 
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e31</penalty>
        <surface>Left</surface>
        <value>10.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
