<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Testlab: 2D L-shape with Dirichlet Boundaries (BL-CSR Format built for GPU, PETSc Linear System Solver)</title>
    <timeloop>TestlabLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>1</output-period>
    <save-final-time>false</save-final-time>
    <format name="VtkHdfV2PostProcessor" />
    <output>
      <variable>U</variable>
    </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>L-shape.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>-5.5</f>
    <solution-comparison-file>poisson_test_ref_L-shape_2D.txt</solution-comparison-file>
    <blcsr>true</blcsr>
    <legacy>false</legacy>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.5</value>
    </dirichlet-boundary-condition>
    <linear-system name="AlienLinearSystem">
      <linear-solver name="PETScSolver">
        <solver name="CG">
          <num-iterations-max>1000</num-iterations-max>
          <stop-criteria-value>1e-8</stop-criteria-value>
          <preconditioner name="Jacobi">
          </preconditioner>
        </solver>
        <verbose>high</verbose>
      </linear-solver>
    </linear-system>
  </fem>
</case>
