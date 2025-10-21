<?xml version="1.0"?>
<case codename="Soildynamics" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Transient traction case</title>
    <timeloop>SoildynamicsLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
     <variable>V</variable>
     <variable>A</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/semi-circle-soil.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>0.08</tmax>
    <dt>0.01</dt>
    <E>6.62e6</E>
    <nu>0.45</nu>
    <rho>2500.0</rho>
    <f>NULL 315.9</f>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <penalty>1.e30</penalty>
    <paraxial-boundary-condition>
      <surface>lower</surface>
    </paraxial-boundary-condition>
    <boundary-conditions>
      <traction>
        <surface>input</surface>
        <traction-input-file>data/semi-circle-soil-traction.txt</traction-input-file>
      </traction>
    </boundary-conditions>
    <result-file>check/test_2D_transient_traction.txt</result-file>
    <linear-system>
      <solver-backend>petsc</solver-backend>
      <preconditioner>ilu</preconditioner>
    </linear-system>
  </fem>
</case>
