<?xml version="1.0"?>
<case codename="Elastoplasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>2D quarter cylinder with two material domains</title>
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
      <filename>meshes/quater_cylinder-multi_domain.msh</filename>
    </mesh>
  </meshes>

  <fem>
    <tmax>21.</tmax>
    <dt>1.</dt>
    <matrix-format>DOK</matrix-format>
    <gp-material-tensor-strategy>global</gp-material-tensor-strategy>

    <material-domain>
      <volume>domain1</volume>
    </material-domain>
    <material-domain>
      <volume>domain2</volume>
    </material-domain>

    <constitutive-law>
      <law>VonMises</law>
      <von-mises>
        <E>70.0e3</E>
        <nu>0.3</nu>
        <sig0>250.</sig0>
      </von-mises>
    </constitutive-law>
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
  </fem>
</case>
