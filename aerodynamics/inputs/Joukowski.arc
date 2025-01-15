<?xml version="1.0"?>
<case codename="aerodynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>aerodynamicsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>u</variable>
     <variable>psi</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/NACA0012.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <farfield-boundary-condition>
      <surface>FarField</surface>
      <angle>0.1</angle>
    </farfield-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>upperAirfoil</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>lowerAirfoil</surface>
      <value>0.0</value>
    </dirichlet-boundary-condition>
  </fem>
</case>
