<?xml version="1.0"?>
<!--
  Case configuration for a Fourier analysis simulation.
  The XML file includes sections for:
    - General simulation settings
    - Mesh configuration details
    - Finite Element Method (FEM) configurations
    - Post-processing options
    - External function definitions
-->
<case codename="Fourier" xml:lang="en" codeversion="1.0">

  <!--
    Arcane-specific settings:
      - title: A descriptive name for the case.
      - timeloop: Defines the specific time-stepping loop used for this Fourier simulation.
  -->
  <arcane>
    <title>Fouriers equation FEM code with maufactured solution</title>
    <timeloop>FourierLoop</timeloop>
  </arcane>

  <!--
    Mesh configuration:
      - filename: Path to the mesh file used in the simulation.
  -->
  <meshes>
    <mesh>
      <filename>meshes/square_-2pi_to_2pi.msh</filename>
    </mesh>
  </meshes>

  <!--
    Function definitions for external assemblies:
      - external-assembly: Specifies the external dynamic link library (DLL) and the class containing the functions.
      - assembly-name: The name of the DLL file containing external functions.
      - class-name: The specific class within the DLL that contains the necessary case functions.
  -->
  <functions>
    <external-assembly>
      <assembly-name>ExternalFunctions.dll</assembly-name>
      <class-name>FemModuleFourier.CaseFunctions</class-name>
    </external-assembly>
  </functions>

  <!--
    FEM (Finite Element Method) settings:
      - boundary-conditions: Specifies boundary conditions including a manufactured solution.
      - manufactured-solution: Configuration for the manufactured solution approach, useful for verification.
        - manufactured-dirichlet: Specifies a Dirichlet boundary condition function.
        - manufactured-source: Specifies a source term function.
        - enforce-Dirichlet-method: Method used to enforce Dirichlet boundary conditions (e.g., Penalty method).
      - lambda: Parameter for the FEM calculations  representing thermal conductivity.
  -->
  <fem>
    <boundary-conditions>
      <manufactured-solution>
       <manufactured-dirichlet function="manufacturedDirichlet">true</manufactured-dirichlet>
       <manufactured-source function="manufacturedSource">true</manufactured-source>
       <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
      </manufactured-solution>
    </boundary-conditions>
    <lambda>1.0</lambda>
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
     <variable>UExact</variable>
   </output>
  </arcane-post-processing>

</case>
