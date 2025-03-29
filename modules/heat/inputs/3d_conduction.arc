<?xml version="1.0"?>
<case codename="Heat" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>3D heat conduction poroblem</title>
    <timeloop>HeatLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>2</output-period>
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/truncated_cube.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>1.75</lambda>
    <tmax>1.</tmax>
    <dt>0.1</dt>
    <Tinit>30.0</Tinit>
    <result-file>check/3d_conduction.txt</result-file>
    <boundary-conditions>
      <dirichlet>
        <surface>top</surface>
        <value>100.6</value>
      </dirichlet>
      <dirichlet>
        <surface>bottom</surface>
        <value>1.6</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>
