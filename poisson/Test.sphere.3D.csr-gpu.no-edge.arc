<?xml version="1.0"?>
<case codename="Poisson" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sphere 3D with CSR sparse matrix format Gpu compatible without edge. The result of this test is compared with poisson_test_ref_sphere_3D.txt</title>
    <timeloop>PoissonLoop</timeloop>
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
      <filename>sphere_cut.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <result-file>poisson_test_ref_sphere_3D.txt</result-file>
    <f>-0.035</f>
    <mesh-type>TETRA4</mesh-type>
    <enforce-Dirichlet-method>Penalty</enforce-Dirichlet-method>
    <penalty>1.e31</penalty>
    <dirichlet-boundary-condition>
      <surface>sphere</surface>
      <value>-5.5</value>
    </dirichlet-boundary-condition>
    <dirichlet-boundary-condition>
      <surface>Cut</surface>
      <value>50.0</value>
    </dirichlet-boundary-condition>
    <linear-system name="HypreLinearSystem">
      <rtol>0.</rtol>
      <atol>1e-5</atol>
      <amg-threshold>0.55</amg-threshold>
    </linear-system>
    <csr-gpu>true</csr-gpu>
    <create-edges>false</create-edges>
  </fem>
</case>