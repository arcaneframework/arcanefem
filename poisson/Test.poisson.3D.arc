<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PoissonLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>L-shape-3D.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>test_3D_L-shape.txt</result-file>
    <mesh-type>TETRA4</mesh-type>
    <dirichlet-boundary-condition>
      <surface>bot</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>bc</surface>
      <value>10.0</value>
    </dirichlet-boundary-condition>
  </fem>
</case>
