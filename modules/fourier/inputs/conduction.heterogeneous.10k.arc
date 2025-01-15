<?xml version="1.0"?>
<!--
  Case configuration for a Fourier analysis simulation.
  This configuration file contains settings for:
    - General simulation parameters
    - Mesh configuration details
    - Finite Element Method (FEM) configurations
    - Post-processing options
-->
<case codename="Fourier" xml:lang="en" codeversion="1.0">

  <!--
    Arcane-specific settings:
      - title: Descriptive name for the simulation case.
      - timeloop: Specifies the time-stepping loop used in this Fourier simulation.
  -->
  <arcane>
    <title>Fouriers equation FEM code with heterogenous material and fine mesh</title>
    <timeloop>FourierLoop</timeloop>
  </arcane>

  <!--
    Mesh configuration:
      - filename: The path to the mesh file used in the simulation.
      - The mesh file is designed for multi-material analysis with a 10,000-element mesh.
  -->
  <meshes>
    <mesh>
      <filename>meshes/multi-material.10k.msh</filename>
    </mesh>
  </meshes>

  <!--
    FEM (Finite Element Method) settings:
      - lambda: Default thermal conductivity or diffusivity coefficient for the simulation.
      - qdot: Specifies the heat source term or volumetric heat generation.
      - boundary-conditions: Defines the boundary conditions for the simulation, including Dirichlet and Neumann conditions.
        - dirichlet: Applies fixed value boundary conditions to specified surfaces.
        - neumann: Applies flux or gradient boundary conditions to specified surfaces.
      - material-property: Specifies the material properties for different volumes in the mesh, allowing for multi-material simulations.
        - volume: Identifies the material volume within the mesh.
        - lambda: Specifies the thermal conductivity or diffusivity for the corresponding material.
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
      - output-period: Determines how often (in simulation steps) the output is generated.
      - format: Specifies the output format, here using VtkHdfV2.
      - output: Defines the variables to be included in the post-processing output.
  -->
  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

</case>
