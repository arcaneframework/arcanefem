<?xml version="1.0"?>
<case codename="Testlab" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Testlab: 2D L-shape with Dirichlet Boundaries (BL-CSR Format built for GPU, PETSc Linear System Solver)</title>
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
      <filename>L-shape-3D.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <f>-5.5</f>
    <blcsr>true</blcsr>
    <legacy>false</legacy>
    <dirichlet-boundary-condition>
      <surface>bot</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>bc</surface>
      <value>10.0</value>
    </dirichlet-boundary-condition>
    <linear-system name="PETScLinearSystem">
      <rtol>0.</rtol>
      <atol>1e-5</atol>
    </linear-system>
  </fem>
</case>

