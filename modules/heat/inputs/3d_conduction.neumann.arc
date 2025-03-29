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
    <boundary-conditions>
      <dirichlet>
        <surface>top</surface>
        <value>100.6</value>
      </dirichlet>
      <neumann>
        <surface>left</surface>
        <value>9.6</value>
      </neumann>
    </boundary-conditions>
  </fem>
</case>
