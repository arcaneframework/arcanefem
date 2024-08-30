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
    <title>Fouriers equation FEM code with quad mesh</title>
    <timeloop>FourierLoop</timeloop>
  </arcane>

  <!--
    Mesh configuration:
      - filename: Path to the mesh file used in the simulation.
  -->
  <meshes>
    <mesh>
      <filename>plancher.quad4.msh</filename>
    </mesh>
  </meshes>

  <!--
    FEM (Finite Element Method) settings:
      - lambda: Thermal conductivity or diffusivity coefficient.
      - qdot: Heat source term or volumetric heat generation.
      - mesh-type: Specifies the type of mesh used in the simulation (e.g., QUAD4).
      - boundary-conditions: Defines the boundary conditions for the simulation.
        - dirichlet: Fixed value boundary condition for specific surfaces.
        - neumann: Flux or gradient boundary condition for specific surfaces.
  -->
  <fem>
    <lambda>1.75</lambda>
    <qdot>1e5</qdot>
    <mesh-type>QUAD4</mesh-type>
    <boundary-conditions>
      <dirichlet>
        <surface>Cercle</surface>
        <value>50.0</value>
      </dirichlet>
      <dirichlet>
        <surface>Bas</surface>
        <value>5.0</value>
      </dirichlet>
      <dirichlet>
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