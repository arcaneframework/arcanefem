<?xml version="1.0"?>
<!--
  Case configuration for a Fourier analysis simulation.
  The XML file includes sections for:
    - General simulation settings
    - Mesh configuration details
    - Finite Element Method (FEM) configurations
    - Post-processing options
-->
<case codename="Fourier" xml:lang="en" codeversion="1.0">

  <!--
    Arcane-specific settings:
      - title: A descriptive name for the case.
      - timeloop: Defines the specific time-stepping loop used for this Fourier simulation.
  -->
  <arcane>
    <title>Fouriers equation FEM code with heterogenous material</title>
    <timeloop>FourierLoop</timeloop>
  </arcane>

  <!--
    Mesh configuration:
      - filename: Path to the mesh file used in the simulation.
      - The mesh file is configured for a multi-material simulation.
  -->
  <meshes>
    <mesh>
      <filename>meshes/multi-material.msh</filename>
    </mesh>
  </meshes>

  <!--
    FEM (Finite Element Method) settings:
      - lambda: Default thermal conductivity or diffusivity coefficient.
      - qdot: Heat source term or volumetric heat generation.
      - boundary-conditions: Defines the boundary conditions for the simulation.
        - dirichlet: Fixed value boundary condition for specific surfaces.
        - neumann: Flux or gradient boundary condition for specific surfaces.
      - material-property: Specifies material properties, like thermal conductivity, for different volumes in the mesh.
        - volume: Name of the material volume within the mesh.
        - lambda: Thermal conductivity or diffusivity for the specified material.
  -->
  <fem>
    <lambda>0.0</lambda>
    <qdot>15.</qdot>
    <boundary-conditions>
      <dirichlet>
        <surface>Left</surface>
        <value>50.0</value>
      </dirichlet>
      <dirichlet>
        <surface>Right</surface>
        <value>5.0</value>
      </dirichlet>
      <neumann>
        <surface>Top</surface>
        <value>0.0</value>
      </neumann>
      <neumann>
        <surface>Bot</surface>
        <value>0.0</value>
      </neumann>
    </boundary-conditions>
    <material-property>
      <volume>Mat1</volume>
      <lambda>100.0</lambda>
    </material-property>
    <material-property>
      <volume>Mat2</volume>
      <lambda>1.0</lambda>
    </material-property>
  </fem>

  <!--
    Post-processing settings:
      - output-period: Defines the frequency (in simulation steps) at which output is generated.
      - format: Specifies the post-processing format, in this case, VtkHdfV2.
      - output: Lists the variables to be output during post-processing.
  -->
  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

</case>
