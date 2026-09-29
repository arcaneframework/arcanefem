<?xml version="1.0"?>
<case codename="Elastoplasticity2" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>2D quarter-cylinder driven by a compiled MFront behaviour (generic law)</title>
    <timeloop>Elastoplasticity2Loop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>1</output-period>
    <output>
      <variable>U</variable>
    </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/quater_cylinder.msh</filename>
    </mesh>
  </meshes>

  <elastoplasticity2>
    <tmax>21.</tmax>
    <dt>1.</dt>
    <constitutive-law>
      <law>MFront</law>
      <mfront>
        <!-- Path to the compiled MFront behaviour library (.so), relative to
             the run directory. Provide the library you compiled with MFront. -->
        <behaviour-file>data/libElastic2D.so</behaviour-file>
        <!-- Symbol prefix / registration name of the behaviour in the library. -->
        <behaviour-name>Elastic2D</behaviour-name>
        <!-- PlaneStrain | PlaneStress | Axisymmetrical | GeneralisedPlaneStrain | Tridimensional -->
        <hypothesis>PlaneStrain</hypothesis>
        <!-- Material properties, in the order declared by the behaviour
             ([YoungModulus, PoissonRatio] here). -->
        <material-properties>70000.0 0.3</material-properties>
      </mfront>
    </constitutive-law>
    <gp-material-tensor-strategy>global</gp-material-tensor-strategy>
    <f>NULL NULL</f>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 NULL</value>
      </dirichlet>
      <dirichlet>
        <surface>bottom</surface>
        <value>NULL 0.0</value>
      </dirichlet>
      <traction>
        <surface>inner</surface>
        <traction-input-file>data/traction_quater_cylinder_20steps.txt</traction-input-file>
      </traction>
    </boundary-conditions>
  </elastoplasticity2>
</case>