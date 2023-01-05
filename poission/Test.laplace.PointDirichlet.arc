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
      <filename>plancher.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>1.75</lambda>
    <qdot>0.0</qdot>
    <result-file>test3_results.txt</result-file>
    <dirichlet-point-condition>
      <node>topLeftCorner</node>
      <value>50.0</value>
    </dirichlet-point-condition>
    <dirichlet-point-condition>
      <node>topRightCorner</node>
      <value>20.0</value>
    </dirichlet-point-condition>
    <dirichlet-point-condition>
      <node>botLeftCorner</node>
      <value>20.0</value>
    </dirichlet-point-condition>
    <dirichlet-point-condition>
      <node>botRightCorner</node>
      <value>50.0</value>
    </dirichlet-point-condition>
    <linear-system name="SequentialBasicLinearSystem" />
  </fem>
</case>
