<?xml version='1.0'?>
<case codename="Passmo" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>Sample</title>
    <timeloop>PassmoLoop</timeloop>
  </arcane>
  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>Displacement</variable>
   </output>
  </arcane-post-processing>

    <mesh>
      <filename>bar.msh</filename>
    </mesh>

    <elastodynamic>
      <max-frequency>5.</max-frequency>

      <input-motion>
        <node-group>1 4 35 36</node-group>
        <acceleration-input-file>sinus1_acc.txt</acceleration-input-file>
        <velocity-input-file>sinus1_vel.txt</velocity-input-file>
        <displacement-input-file>sinus1_disp.txt</displacement-input-file>
        <rigid-base>true</rigid-base>true>
        <x-component>true</x-component>
        <y-component>false</y-component>
        <z-component>false</z-component>
        <amplification-factors>1. 0. 0.</amplification-factors>
      </input-motion>
      <analysis-type>planestrain</analysis-type>
      <start>0.</start>
      <final-time>10.</final-time>
    <simple name="deltat" type="real">
      <description>Timestep value for simulation</description>
    </simple>
    <simple name="beta" type="real" default="0.25" optional="true">
      <description>Newmark Beta coefficient</description>
    </simple>
    <simple name="gamma" type="real" default="0.5" optional="true">
      <description>Newmark Gamma coefficient</description>
    </simple>
    <simple name="alfam" type="real" default="0." optional="true">
      <description>Coefficient related to mass terms in Newmark Generalized alfa-method</description>
    </simple>
    <simple name="alfaf" type="real" default="0." optional="true">
      <description>Coefficient related to force terms in Newmark Generalized alfa-method</description>
    </simple>
    <simple name="is_alfa_method" type="bool" default="false" optional="true">
      <description>Boolean which is true if Newmark Generalized alfa-method is used</description>
    </simple>
    <simple name="nb_dofs_per_node" type="integer" default="3" optional="true">
      <description>Number of dofs per node</description>
    </simple>    </elastodynamic>
</case>
