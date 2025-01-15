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
    <title>Fouriers equation FEM code fine mesh</title>
    <timeloop>FourierLoop</timeloop>
  </arcane>

  <!--
    Mesh configuration:
      - filename: The path to the mesh file used in the simulation.
  -->
  <meshes>
    <mesh>
      <filename>meshes/plancher.10k.msh</filename>
    </mesh>
  </meshes>

  <!--
    FEM (Finite Element Method) settings:
      - lambda: Thermal conductivity or diffusivity coefficient.
      - qdot: Heat source term or volumetric heat generation.
      - boundary-conditions: Defines the boundary conditions for the simulation.
        - dirichlet: Fixed value boundary condition with penalty enforcement for specified surfaces.
        - neumann: Flux or gradient boundary condition for specified surfaces.
  -->
  <fem>
    <lambda>1.75</lambda>
    <qdot>1e5</qdot>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e12</penalty>
        <surface>Cercle</surface>
        <value>50.0</value>
      </dirichlet>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e12</penalty>
        <surface>Bas</surface>
        <value>5.0</value>
      </dirichlet>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e12</penalty>
        <surface>Haut</surface>
        <value>21.0</value>
      </dirichlet>
      <neumann>
        <surface>Droite</surface>
        <value>15.0</value>
      </neumann>
      <neumann>
        <surface>Gauche</surface>
        <value>0.0</value>
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
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

</case>
