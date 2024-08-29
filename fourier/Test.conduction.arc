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
      <filename>plancher.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>1.75</lambda>
    <qdot>1e5</qdot>
    <result-file>test1_results.txt</result-file>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e12</penalty>
        <surface>Cercle</surface>
        <value>50.0</value>
      </dirichlet>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e12</penalty>
        <surface>Bas</surface>
        <value>5.0</value>
      </dirichlet>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e12</penalty>
        <surface>Haut</surface>
        <value>21.0</value>
      </dirichlet>
      <neumann>
        <surface>Droite</surface>
        <value>15.0</value>
      </neumann>
      <neumann>
        <surface>Gauche</surface>
        <value>0.0</value>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
