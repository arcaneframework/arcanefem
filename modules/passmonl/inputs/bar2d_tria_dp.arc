<?xml version='1.0'?>
<case codename="Passmonl" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PassmonlLoop</timeloop>
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

  <n-l-dynamic>
    <analysis-type>planestrain</analysis-type>
    <start>0.</start>
    <final-time>2.0</final-time>
    <deltat>0.08</deltat>

    <init-elast-type>lame</init-elast-type>
    <Itemax>30</Itemax>

      <boundary-conditions>
          <dirichlet>
              <surface>surfaceleft</surface>
              <value>0.0 0.0 0.0</value>
          </dirichlet>

          <dirichlet>
              <surface>surfaceright</surface>
              <value>1.0 0.0 0.0</value>
          </dirichlet>
      </boundary-conditions>

    <nonlin-algo-type>modnewtonraphson</nonlin-algo-type>
    <integration-type>femcell</integration-type>
    <law-input-param>dp_params.txt</law-input-param>
    <law-model>
       <cell-group>volume</cell-group>
       <law-type>Druckerprager</law-type>
       <nb-law-param>7</nb-law-param>
       <nb-law-hist-param>1</nb-law-hist-param>
       <i-law-param>0</i-law-param>
    </law-model>

  </n-l-dynamic>
</case>
