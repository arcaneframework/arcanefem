<?xml version="1.0"?>
<case codename="Elastoplasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>2D Cartesian cantilever beam: von Mises loading and unloading</title>
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
      <generator name="Cartesian2D">
        <nb-part-x>0</nb-part-x>
        <nb-part-y>1</nb-part-y>
        <origin>0.0 0.0</origin>
        <generate-sod-groups>true</generate-sod-groups>
        <x><n>40</n><length>5.0</length></x>
        <y><n>8</n><length>1.0</length></y>
      </generator>
    </mesh>
  </meshes>

  <fem>
    <!-- The module solves t=1,...,20 for tmax=21 and dt=1. -->
    <tmax>21.</tmax>
    <dt>1.</dt>
    <hex-quad-mesh>true</hex-quad-mesh>
    <matrix-format>BSR</matrix-format>
    <gp-material-tensor-strategy>global</gp-material-tensor-strategy>
    <constitutive-law>
      <law>VonMises</law>
      <von-mises>
        <E>70.0e3</E>
        <nu>0.3</nu>
        <sig0>250.</sig0>
      </von-mises>
    </constitutive-law>
    <linear-system name="HypreLinearSystem">
      <solver>cg</solver>
      <atol>1e-10</atol>
      <rtol>0</rtol>
    </linear-system>
    <f>NULL NULL</f>
    <boundary-conditions>
      <dirichlet>
        <surface>XMIN</surface>
        <value>0.0 0.0</value>
      </dirichlet>
      <traction>
        <surface>YMAX</surface>
        <!-- Positive table pressure acts inward, hence downward on YMAX. -->
        <traction-input-file>data/pressure_beam_20steps.txt</traction-input-file>
      </traction>
    </boundary-conditions>
  </fem>
</case>
