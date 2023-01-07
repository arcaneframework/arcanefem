<?xml version="1.0"?>
<case codename="Elastodynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ElastodynamicsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>10</output-period>
   <output>
     <variable>U</variable>
     <variable>V</variable>
     <variable>A</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>bar.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>0.1</tmax>
    <dt>0.01</dt>
    <etam>0.01</etam>
    <etak>0.01</etak>
    <alpm>0.20</alpm>
    <alpf>0.40</alpf>
    <rho>1.0</rho>
    <E>1000.0</E>
    <nu>0.3</nu>
    <f2>-1.0</f2>
    <dirichlet-boundary-condition>
      <surface>left</surface>
      <u1>0.0</u1>
      <u2>0.0</u2>
    </dirichlet-boundary-condition>
  </fem>
</case>
