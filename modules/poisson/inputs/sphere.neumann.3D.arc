<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sphere 3D</title>
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
      <filename>meshes/sphere_cut.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>check/poisson_test_ref_sphere_neumann_3D.txt</result-file>
    <f>5.5</f>
    <boundary-conditions>
      <dirichlet>
        <surface>horizontal</surface>
        <value>0.5</value>
      </dirichlet>
      <neumann>
        <surface>curved</surface>
        <value>0.35 1.65 3.75</value>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
