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
     <variable>Flux</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>plate.fine.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>1.75</lambda>
    <tmax>20.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
    <dirichlet-boundary-condition>
      <surface>left</surface>
      <value>10.0</value>
    </dirichlet-boundary-condition>
  </fem>
</case>
