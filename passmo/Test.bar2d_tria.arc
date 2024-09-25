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
      <filename>bar_dynamic.msh</filename>
      <initialization>
        <variable><name>Rho</name><value>1.000000</value><group>volume</group></variable>
        <variable><name>Lambda</name><value>576.9230769</value><group>volume</group></variable>
        <variable><name>Mu</name><value>384.6153846</value><group>volume</group></variable>
      </initialization>
    </mesh>
  </meshes>

  <elastodynamic>
    <analysis-type>planestrain</analysis-type>
    <start>0.</start>
    <final-time>2.0</final-time>
    <deltat>0.08</deltat>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>

    <init-elast-type>lame</init-elast-type>

    <dirichlet-surface-condition>
      <surface>surfaceleft</surface>
      <Ux>0.0</Ux>
      <Uy>0.0</Uy>
    </dirichlet-surface-condition>

    <dirichlet-surface-condition>
      <surface>surfaceright</surface>
      <Ux>1.0</Ux>
    </dirichlet-surface-condition>

  </elastodynamic>
</case>
