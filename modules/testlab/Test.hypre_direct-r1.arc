<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Testlab: 2D L-shape with Dirichlet Boundaries (BL-CSR Format built for GPU, Hypre Linear System Solver)</title>
    <timeloop>TestlabLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>1</output-period>
    <save-final-time>false</save-final-time>
    <output>
      <variable>U</variable>
    </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>L-shape-r1.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>-5.5</f>
    <blcsr>true</blcsr>
    <legacy>false</legacy>
    <dirichlet-boundary-condition>
      <surface>boundary</surface>
      <value>0.5</value>
    </dirichlet-boundary-condition>
      <rtol>0.</rtol>
      <atol>1e-5</atol>
      <amg-threshold>0.25</amg-threshold>
    </linear-system>
  </fem>
</case>
