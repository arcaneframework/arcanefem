<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Testlab: 2D Porous medium with  multiple Dirichlet Boundary (DOK Format, PETSc Linear System Solver)</title>
    <timeloop>TestlabLoop</timeloop>
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
      <filename>porous-medium.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>-5.5</f>
    <dirichlet-boundary-condition>
      <surface>wall</surface>
      <value>0.5</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore1</surface>
      <value>300.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore2</surface>
      <value>250.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore3</surface>
      <value>230.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore4</surface>
      <value>220.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore5</surface>
      <value>283.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore6</surface>
      <value>273.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore7</surface>
      <value>265.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore8</surface>
      <value>241.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore9</surface>
      <value>224.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore10</surface>
      <value>256.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>pore11</surface>
      <value>267.0</value>
    </dirichlet-boundary-condition>
    <linear-system>
      <solver-backend>petsc</solver-backend>
    </linear-system>
  </fem>
</case>
