<?xml version='1.0'?>
<case codename="Passmo" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PassmoLoop</timeloop>
  </arcane>
  <arcane-post-processing>
<arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>Displ</variable>
     <variable>Vel</variable>
     <variable>Acc</variable>
   </output>
   <!-- <ensight7gold>
    <binary-file>false</binary-file>
   </ensight7gold>-->
 </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>barq4dbg.msh</filename>
      <initialization>
        <variable><name>Rho</name><value>2200.0</value><group>surface</group></variable>
        <variable><name>Young</name><value>6.e7</value><group>surface</group></variable>
        <variable><name>Nu</name><value>0.3</value><group>surface</group></variable>
      </initialization>
    </mesh>
  </meshes>

  <elastodynamic>
    <analysis-type>planestrain</analysis-type>
    <start>0.</start>
    <final-time>0.05</final-time>
    <deltat>0.01</deltat>
    <beta>0.25</beta>
    <gamma>0.5</gamma>
    <alfa_method>false</alfa_method>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <penalty>1.e30</penalty>
    <linop-nstep>10</linop-nstep>
    <gz>10.0</gz>

    <init-elast-type>young</init-elast-type>

    <dirichlet-boundary-condition>
      <surface>left</surface>
      <Ux>0.0</Ux>
      <Uy>0.0</Uy>
    </dirichlet-boundary-condition>

    <neumann-boundary-condition>
      <surface>right</surface>
      <Y-val>0.01</Y-val>
    </neumann-boundary-condition>

    <linear-system name="SequentialBasicLinearSystem" />
  </elastodynamic>
</case>
