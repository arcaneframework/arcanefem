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
      <filename>meshes/plate.quad.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <lambda>1.75</lambda>
    <tmax>20.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
    <result-file>check/2d_conduction_neumann_pointBC.quad.txt</result-file>
    <boundary-conditions>
      <dirichlet-point>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e31</penalty>
        <node>topLeft</node>
        <value>1.8</value>
      </dirichlet-point>
      <dirichlet-point>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e31</penalty>
        <node>botRight</node>
        <value>31.0</value>
      </dirichlet-point>
      <neumann>
        <surface>left</surface>
        <value>9.6</value>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
