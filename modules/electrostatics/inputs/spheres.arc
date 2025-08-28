<?xml version="1.0"?>
<case codename="Electrostatics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ElectrostaticsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>Phi</variable>
     <variable>E</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/spheres_in_sphere.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <rho>0.0</rho>
    <epsilon>1.0</epsilon>

    <boundary-conditions>
      <dirichlet>
        <surface>outerSphere3</surface>
        <value>-1.0</value>
      </dirichlet>
      <dirichlet>
        <surface>outerSphere2</surface>
        <value>-1.0</value>
      </dirichlet>
      <dirichlet>
        <surface>outerSphere1</surface>
        <value>-1.0</value>
      </dirichlet>
      <dirichlet>
        <surface>outerSphere</surface>
        <value>0.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
