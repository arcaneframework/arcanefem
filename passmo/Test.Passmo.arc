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
      <filename>bar_dynamic.msh</filename>
    <init>
      <variable name="Density" value="1." group="volume" />
      <variable nom="YoungModulus" value="591.715976" group="volume" />
      <variable nom="PoissonRatio" value="0.3" group="volume" />
    </init>
    </mesh>

    <elastodynamic>
      <analysis-type>planestrain</analysis-type>
      <start>0.</start>
      <final-time>3.</final-time>
      <deltat>0.08</detltat>
      <beta>0.25</beta>
      <gamma>0.5</gamma>
      <is_alfa_method>false</is_alfa_method>
      <nb_dofs_per_node>2</nb_dofs_per_node>

      <dirichlet-boundary-condition>
        <surface>surfaceleft</surface>
        <type>DisplacementX</type>
        <constant-value>0.0</constant-value>
      </dirichlet-boundary-condition>

      <dirichlet-boundary-condition>
        <surface>surfaceleft</surface>
        <type>DisplacementY</type>
        <constant-value>0.0</constant-value>
      </dirichlet-boundary-condition>

    <neumann-boundary-condition>
      <surface>surfaceright</surface>
      <X-constant>0.0</X-constant>
      <Y-file>input_traction.txt</Y-file>
    </neumann-boundary-condition>

    <linear-system>
      <solver-backend>petsc</solver-backend>
      <preconditioner>ilu</preconditioner>
    </linear-system>
  </elastodynamic>
</case>
