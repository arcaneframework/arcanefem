<?xml version="1.0"?>
<case codename="Fourier" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>FourierLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
     <variable>UExact</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>square_-2pi_to_2pi.msh</filename>
    </mesh>
  </meshes>

  <functions>
    <external-assembly>
      <assembly-name>ExternalFunctions.dll</assembly-name>
      <class-name>FemModule.CaseFunctions</class-name>
    </external-assembly>
  </functions>

  <fem>
    <manufactured-dirichlet-condition function="manufacturedDirichlet">true</manufactured-dirichlet-condition>
    <manufactured-source-condition function="manufacturedSource">true</manufactured-source-condition>
    <lambda>1.0</lambda>
  </fem>
</case>
