<?xml version="1.0"?>
<!--
  Case configuration for an acoustics simulation. 
  The following sections define:
    - General simulation settings
    - Mesh configurations details
    - Finite Element Method (FEM) configurations
    - Post-processing options
-->
<case codename="Acoustics" xml:lang="en" codeversion="1.0">

  <!--
    Arcane-specific settings:
      - title: The name or description of the case.
      - timeloop: Defines the time-stepping loop for the simulation.
  -->
  <arcane>
    <title>Submarine Toy Case</title>
    <timeloop>AcousticsLoop</timeloop>
  </arcane>


  <!--
    Mesh configurations:
      - filename: Path to the mesh file used in the simulation.
  -->
  <meshes>
    <mesh>
      <filename>meshes/sub.quad.msh</filename>
    </mesh>
  </meshes>


  <!--
    FEM (Finite Element Method) settings:
      - hex-quad-mesh: Indicates if the mesh is a hexahedral or quadrilateral mesh.
      - kc2: Coefficient used in the FEM calculations.
      - boundary-conditions: Defines neumann boundary conditions for the simulation.
      - linear-system: Specifies the linear system solver to use.
      - result-file: File for validation (optional)
  -->
  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <kc2>.11e1</kc2>
    <boundary-conditions>
      <neumann>
        <surface>inner1</surface>
        <value>1.0</value>
      </neumann>
    </boundary-conditions>
    <linear-system name="SequentialBasicLinearSystem" />
    <result-file>check/sub_2D.quad.txt</result-file>
  </fem>

  <!--
    Post-processing settings:
      - output-period: Defines how often output should be generated.
      - output: Specifies which variables are to be outputted.
  -->
  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

</case>
