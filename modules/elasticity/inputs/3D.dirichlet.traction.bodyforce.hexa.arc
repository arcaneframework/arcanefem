<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>ElasticityLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/truncated_cube.hexa.msh</filename>
      <subdivider>
        <nb-subdivision>0</nb-subdivision>
      </subdivider>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <result-file>check/3D_dirichlet_traction_bodyforce.hexa.txt</result-file>
    <E>200e9</E>
    <nu>0.3</nu>
    <f>-9.8e1, -7.5e1, 5.9e1</f>
    <boundary-conditions>
      <dirichlet>
        <surface>top</surface>
        <value>0. 0. 0.</value>
      </dirichlet>
      <traction>
        <surface>bottom</surface>
        <value>12.9e11, -14.5e11, -18.8e11</value>
      </traction>
    </boundary-conditions>
  </fem>
</case>