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
<case codename="FourierNL" xml:lang="en" codeversion="1.0">

  <!--
    Arcane-specific settings:
      - title: A descriptive name for the case.
      - timeloop: Defines the specific time-stepping loop used for this Fourier simulation.
  -->
  <arcane>
    <title>Fouriers equation FEM code with maufactured solution</title>
    <timeloop>FourierNLLoop</timeloop>
  </arcane>

  <!--
    Mesh configuration:
      - filename: Path to the mesh file used in the simulation.
  -->
  <meshes>
      <mesh>
        <generator name="Cartesian3D" >
          <nb-part-x>1</nb-part-x>
          <nb-part-y>1</nb-part-y>
          <nb-part-z>1</nb-part-z>
          <origin>0.0 0.0 0.0</origin>
          <generate-sod-groups>true</generate-sod-groups>
          <x><n>10</n><length>1.0</length></x>
          <y><n>4</n><length>1.0</length></y>
          <z><n>4</n><length>1.0</length></z>
        </generator>
      </mesh>
    </meshes>

  <!--
    FEM (Finite Element Method) settings:
      - boundary-conditions: Specifies boundary conditions
  -->
  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <boundary-conditions>
        <dirichlet>
            <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
            <surface>XMIN</surface>
            <value>0.0</value>
        </dirichlet>
        <dirichlet>
            <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
            <surface>XMAX</surface>
            <value>1.0</value>
        </dirichlet>
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
   <output>
     <variable>U</variable>
     <variable>UExact</variable>
   </output>
  </arcane-post-processing>

</case>
