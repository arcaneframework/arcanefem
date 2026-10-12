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
    <timeloop>FourierHybridLoop</timeloop>
  </arcane>

  <!--
    Mesh configuration:
      - filename: Path to the mesh file used in the simulation.
      - The mesh file is configured for a multi-material simulation.
  -->
  

  <!--
  <meshes>
    <mesh>
      <filename>my_cases/meshes/skyscraper.msh</filename>
    </mesh>
  </meshes>
  -->

  <meshes>
    <mesh>
      <generator name="Cartesian2D">
        <x>
          <n>200</n>
          <length>1.0</length>
        </x>

        <y>
          <n>200</n>
          <length>1.0</length>
        </y>

        <origin>0.0 0.0</origin>

        <nb-part-x>0</nb-part-x>
        <nb-part-y>0</nb-part-y>

        <!-- <generate-sod-groups>true</generate-sod-groups> -->

      </generator>
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
    <lambda>1.0</lambda>

    <qdot>15.0</qdot>

    <matrix-format>DOK</matrix-format>
	
    <use-skyscraper>true</use-skyscraper>
    <linear-system name="HypreLinearSystem">
    </linear-system>

    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>     
        <surface>YMIN</surface>
        <value>10.0</value>
      </dirichlet>

      <dirichlet>
        <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
        <surface>YMAX</surface>
        <value>50.0</value>
      </dirichlet>

      <neumann>
        <surface>XMAX</surface>
        <value>0.0</value>
      </neumann>

      <neumann>
        <surface>XMIN</surface>
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
   <!-- <format name="VtkHdfV2PostProcessor" /> -->
   <output>
     <variable>U</variable>
     <variable>CellLambda</variable>
   </output>
  </arcane-post-processing>

</case>
