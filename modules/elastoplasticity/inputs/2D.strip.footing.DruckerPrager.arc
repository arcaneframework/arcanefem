<?xml version="1.0"?>
<case codename="Elastoplasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>2D strip footing geomechanics test from PSD</title>
    <timeloop>ElastoplasticityLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <filename>meshes/strip_footing.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>13.</tmax>
    <dt>1.</dt>
    <constitutive-law>
      <law>DruckerPrager</law>
      <drucker-prager>
        <E>1.0e7</E>
        <nu>0.48</nu>
        <cohesion>450.0</cohesion>
        <friction-angle>0.35</friction-angle>
      </drucker-prager>
    </constitutive-law>
    <gp-material-tensor-strategy>global</gp-material-tensor-strategy>
    <f>NULL NULL</f>
    <boundary-conditions>
      <dirichlet>
        <surface>left</surface>
        <value>0.0 NULL</value>
        <enforce-Dirichlet-method>RowColumnElimination</enforce-Dirichlet-method>
      </dirichlet>
      <dirichlet>
        <surface>right</surface>
        <value>0.0 NULL</value>
      </dirichlet>
      <dirichlet>
        <surface>bottom</surface>
        <value>NULL 0.0</value>
      </dirichlet>
      <dirichlet>
        <surface>footing</surface>
        <value>NULL 0.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>