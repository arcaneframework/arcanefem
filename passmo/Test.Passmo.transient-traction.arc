<?xml version='1.0'?>
<case codename="Passmo" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PassmoLoop</timeloop>
  </arcane>
  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>displ</variable>
     <variable>vel</variable>
     <variable>acc</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>semi-circle-soil.msh</filename>
      <initialization>
        <variable><name>rho</name><value>2500.0</value><group>soil</group></variable>
        <variable><name>young</name><value>6.62e6</value><group>soil</group></variable>
        <variable><name>nu</name><value>0.45</value><group>soil</group></variable>
      </initialization>
    </mesh>
  </meshes>

  <elastodynamic>
    <analysis-type>planestrain</analysis-type>
    <start>0.</start>
    <final-time>0.08</final-time>
    <deltat>0.01</deltat>
    <beta>0.25</beta>
    <gamma>0.5</gamma>
    <alfa_method>false</alfa_method>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <penalty>1.e30</penalty>
    <linop-nstep>10</linop-nstep>
    <gz>10.0</gz>

    <paraxial-boundary-condition>
      <surface>lower</surface>
    </paraxial-boundary-condition>

    <neumann-boundary-condition>
      <surface>input</surface>
      <X-val>0.01</X-val>
      <Y-val>0.01</Y-val>
    </neumann-boundary-condition>

    <linear-system>
      <solver-backend>petsc</solver-backend>
      <preconditioner>ilu</preconditioner>
    </linear-system>
  </elastodynamic>
</case>
