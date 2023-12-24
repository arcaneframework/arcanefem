<?xml version='1.0'?>
<case codename="Passmo" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PassmoLoop</timeloop>
  </arcane>
  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>Displ</variable>
   </output>
 </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>bardbg.msh</filename>
      <initialization>
        <variable><name>Rho</name><value>2200.000000</value><group>surface</group></variable>
        <variable><name>Young</name><value>6.000000e7</value><group>surface</group></variable>
        <variable><name>Nu</name><value>0.300000000000000</value><group>surface</group></variable>
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
    <penalty>1.e64</penalty>
    <linop-nstep>10</linop-nstep>

    <init-elast-type>young</init-elast-type>

    <dirichlet-boundary-condition>
      <surface>left</surface>
      <Ux>0.0</Ux>
      <Uy>0.0</Uy>
    </dirichlet-boundary-condition>
    
    <nint1>1</nint1>
    <nint2>1</nint2>

    <dirichlet-boundary-condition>
      <surface>right</surface>
      <Ux>1.0</Ux>
    </dirichlet-boundary-condition>

    <linear-system name="SequentialBasicLinearSystem" />
  </elastodynamic>
</case>
