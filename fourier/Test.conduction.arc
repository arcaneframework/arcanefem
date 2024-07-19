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
    <enforce-Dirichlet-method>WeakPenalty</enforce-Dirichlet-method>
    <penalty>1.e12</penalty>
    <dirichlet-boundary-condition>
      <surface>Cercle</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Bas</surface>
      <value>5.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Haut</surface>
      <value>21.0</value>
    </dirichlet-boundary-condition>
    <neumann-boundary-condition>
      <surface>Droite</surface>
      <value>15.0</value>
    </neumann-boundary-condition>
    <neumann-boundary-condition>
      <surface>Gauche</surface>
      <value>0.0</value>
    </neumann-boundary-condition>
  </fem>
</case>
