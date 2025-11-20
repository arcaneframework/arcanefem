<?xml version="1.0"?>
<case codename="Elasticity" xml:lang="en" codeversion="1.0">
  <arcane>
    <title>2D Linear Elasticity With Cartesian Mesh</title>
    <timeloop>ElasticityLoop</timeloop>
  </arcane>

  <arcane-post-processing>
   <output-period>1</output-period>
   <output>
     <variable>U</variable>
   </output>
  </arcane-post-processing>

  <meshes>
    <mesh>
      <generator name="Cartesian2D" >
        <nb-part-x>2</nb-part-x> 
        <nb-part-y>1</nb-part-y>
        <origin>0.0 0.0</origin>
        <generate-sod-groups>true</generate-sod-groups>
        <x><n>10</n><length>1.0</length></x>
        <y><n>2</n><length>0.1</length></y>
      </generator>
    </mesh>
  </meshes>

  <fem>
    <hex-quad-mesh>true</hex-quad-mesh>
    <E>21.0e5</E>
    <nu>0.28</nu>
    <f>-1.0 NULL</f>
    <boundary-conditions>
      <dirichlet>
        <surface>XMIN</surface>
        <value>0.0 0.0</value>
      </dirichlet>
      <dirichlet>
        <surface>XMAX</surface>
        <value>NULL 1.0</value>
      </dirichlet>
    </boundary-conditions>
  </fem>
</case>