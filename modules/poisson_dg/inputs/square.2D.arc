<?xml version="1.0"?>
<case codename="Poisson_dg" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Cut circle 2D with DG solver</title>
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
      <filename>meshes/square_tria.med</filename>
    </mesh>
  </meshes>

  <fem>
    <f>5.5</f>
  </fem>
</case>
