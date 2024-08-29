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
    <boundary-conditions>
      <dirichlet>
        <surface>Left</surface>
        <value>50.0</value>
      </dirichlet>
      <dirichlet>
        <surface>Right</surface>
        <value>5.0</value>
      </dirichlet>
      <neumann>
        <surface>Top</surface>
        <value>0.0</value>
      </neumann>
      <neumann>
        <surface>Bot</surface>
        <value>0.0</value>
      </neumann>
    </boundary-conditions>
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
