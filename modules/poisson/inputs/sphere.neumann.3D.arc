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
        <valueX>0.35</valueX>
        <valueY>1.65</valueY>
        <valueZ>3.75</valueZ>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
