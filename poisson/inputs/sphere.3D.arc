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
    <result-file>check/poisson_test_ref_sphere_3D.txt</result-file>
    <mesh-type>TETRA4</mesh-type>
    <boundary-conditions>
      <dirichlet>
        <surface>horizontal</surface>
        <value>0.5</value>
      </dirichlet>
    </boundary-conditions>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <epsilon>1e-15</epsilon>
    </linear-system>
  </fem>
</case>
