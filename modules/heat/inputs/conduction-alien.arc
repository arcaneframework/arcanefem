<?xml version="1.0"?>
<case codename="Heat" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>HeatLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>2</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>NodeTemperature</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/plate.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <lambda>1.75</lambda>
    <tmax>20.</tmax>
    <dt>0.4</dt>
    <Tinit>30.0</Tinit>
    <boundary-conditions>
      <dirichlet>
        <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
        <penalty>1.e31</penalty>
        <surface>left</surface>
        <value>10.0</value>
      </dirichlet>
    </boundary-conditions>
    <linear-system name="AlienLinearSystem">
      <linear-solver name="PETScSolver">
        <exec-space>Device</exec-space>
        <memory-type>Host</memory-type>
        <solver name="BiCGStab">
          <num-iterations-max>1000</num-iterations-max>
          <stop-criteria-value>1e-8</stop-criteria-value>
          <preconditioner name="GAMG">
            <gamg-type>classical</gamg-type>
            <gamg-threshold>0.15</gamg-threshold>
            <gamg-max-levels>25</gamg-max-levels>
            <gamg-agg-nsmooths>1</gamg-agg-nsmooths>
            <gamg-aggressive-coarsening>1</gamg-aggressive-coarsening>
            <gamg-aggressive-square-graph>1</gamg-aggressive-square-graph>
          </preconditioner>
        </solver>
      </linear-solver>
    </linear-system>
  </fem>
</case>

