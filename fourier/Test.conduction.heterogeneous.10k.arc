<?xml version="1.0"?>
<case codename="Fourier" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>FourierLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>multi-material.10k.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>0.0</lambda>
    <qdot>15.</qdot>
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
