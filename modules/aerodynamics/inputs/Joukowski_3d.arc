<?xml version="1.0"?>
<case codename="aerodynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>aerodynamicsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>u</variable>
     <variable>psi</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/aerodynamics_3d_coarse.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <farfield-boundary-condition>
      <surface>outer</surface>
      <angle>0.1</angle>
    </farfield-boundary-condition>
    <boundary-conditions>
      <dirichlet>
        <surface>airplane</surface>
        <value>0.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
