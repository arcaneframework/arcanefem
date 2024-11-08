<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PoissonLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>0</output-period>
   <format name="VtkHdfV2PostProcessor" />
   <output>
     <variable>U</variable>
   </output>
   <save-final-time>false</save-final-time>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>sphere_cut.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>poisson_test_ref_sphere_3D.txt</result-file>
    <f>5.5</f>
    <mesh-type>TETRA4</mesh-type>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <penalty>1.e31</penalty>
    <dirichlet-boundary-condition>
      <surface>horizontal</surface>
      <value>0.5</value>
    </dirichlet-boundary-condition>
    <linear-system name="HypreLinearSystem">
      <rtol>0.</rtol>
      <atol>1e-5</atol>
      <max-iter>10</max-iter>
      <amg-threshold>0.55</amg-threshold>
      <verbosity>4</verbosity>
      <amg-coarsener>8</amg-coarsener>
    </linear-system>
    <csr>true</csr>
    <create-edges>false</create-edges>
  </fem>
</case>
