<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Testlab: 2D Circle with Curved Neumann and Horizontal Dirichlet Boundaries (CSR Format, Hypre Linear System Solver)</title>
    <timeloop>TestlabLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>circle_cut.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <csr>true</csr>
    <solution-comparison-file>poisson_test_ref_circle_trac_2D.txt</solution-comparison-file>
    <f>5.5</f>
    <dirichlet-boundary-condition>
      <surface>horizontal</surface>
      <value>0.5</value>
    </dirichlet-boundary-condition>
    <neumann-boundary-condition>
      <surface>vertical</surface>
      <value>1e4</value>
    </neumann-boundary-condition>
    <linear-system name="HypreLinearSystem">
      <rtol>0.</rtol>
      <atol>1e-5</atol>
      <amg-threshold>0.25</amg-threshold>
    </linear-system>
  </fem>
</case>
