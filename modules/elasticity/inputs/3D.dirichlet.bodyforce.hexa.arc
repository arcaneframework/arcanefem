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
    <result-file>check/3D_dirichlet_bodyforce.hexa.txt</result-file>
    <hex-quad-mesh>true</hex-quad-mesh>
    <E>200e9</E>
    <nu>0.3</nu>
    <f>-9.8e12, -7.5e12, 5.9e12</f>
    <boundary-conditions>
      <dirichlet>
        <surface>top</surface>
        <value>1. 2. 8.</value>
      </dirichlet>
      <dirichlet>
        <surface>bottom</surface>
        <value>12.9, -14.5, -18.8</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>