<?xml version="1.0"?>
<!--
  Case configuration for a Fourier analysis simulation.
  This file includes settings for:
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
    <title>Fouriers equation FEM code</title>
    <timeloop>FourierLoop</timeloop>
  </arcane>

  <!--
    Mesh configuration:
      - filename: The path to the mesh file used in the simulation.
  -->
  <meshes>
    <mesh>
      <filename>meshes/boolean.msh</filename>
    </mesh>
  </meshes>

  <!--
    FEM (Finite Element Method) settings:
      - lambda: Thermal conductivity or diffusivity coefficient.
      - qdot: Heat source term or volumetric heat generation.
      - result-file: File where simulation results will be saved.
      - boundary-conditions: Defines the boundary conditions for the simulation.
        - dirichlet: Fixed value boundary condition with penalty enforcement for specified surfaces.
        - neumann: Flux or gradient boundary condition for specified surfaces.
  -->
  <fem>
    <lambda>23.5</lambda>
    <qdot>1.123e-2</qdot>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <surface>left</surface>
        <value>55.0</value>
      </dirichlet>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <surface>right</surface>
        <value>12.0</value>
      </dirichlet>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <surface>top</surface>
        <value>166.0</value>
      </dirichlet>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <surface>bottom</surface>
        <value>32.0</value>
      </dirichlet>
      <neumann>
        <surface>front</surface>
        <value>56.0</value>
      </neumann>
      <neumann>
        <surface>back</surface>
        <value>3.67</value>
      </neumann>
    </boundary-conditions>
  </fem>

  <!--
    Post-processing settings:
      - output-period: Defines how often (in simulation steps) the output is generated.
      - format: Specifies the post-processing format, here using VtkHdfV2.
      - output: Defines the variables to be included in the post-processing output.
  -->
  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

</case>
