<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sphere 3D</title>
    <timeloop>PoissonLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/sphere_cut.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>5.5</f>
    <boundary-conditions>
      <dirichlet>
        <surface>horizontal</surface>
        <value>0.5</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
