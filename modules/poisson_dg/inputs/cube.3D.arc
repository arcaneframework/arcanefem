<?xml version="1.0"?>
<case codename="Poisson_dg" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Constant patch test on a 3D cube</title>
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
      <filename>meshes/3x3x3_cube_hexa8.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>0.0</f>
    <penalty>12.0</penalty>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>1.0</value>
      </dirichlet>
      <dirichlet>
        <surface>right</surface>
        <value>0.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
