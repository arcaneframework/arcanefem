<?xml version="1.0"?>
<case codename="FemTest" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>FemTestLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>multi-material.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>0.0</lambda>
    <qdot>15.</qdot>
    <result-file>test2_results.txt</result-file>
    <dirichlet-boundary-condition>
      <surface>Left</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Right</surface>
      <value>5.0</value>
    </dirichlet-boundary-condition>
    <neumann-boundary-condition>
      <surface>Top</surface>
      <value>0.0</value>
    </neumann-boundary-condition>
    <neumann-boundary-condition>
      <surface>Bot</surface>
      <value>0.0</value>
    </neumann-boundary-condition>
    <material-property>
      <volume>Mat1</volume>
      <lambda>100.0</lambda>
    </material-property>
    <material-property>
      <volume>Mat2</volume>
      <lambda>1.0</lambda>
    </material-property>
  </fem>
</case>
