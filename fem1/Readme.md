### Notes on the solver ###

The code here is a simple FEM code used to solve a conduction problem on unstructured 2D mesh. 

### Mesh ### 

The mesh `plancher.msh` is provided in the `FemTest1.arc` file 

```xml
  <meshes>
    <mesh>
      <filename>plancher.msh</filename>
    </mesh>
  </meshes>
```

Please not that use version 4 `.msh` file from `Gmsh`. 

### Boundary conditions ###

These are setup via the `Fem1Module::_initBoundaryconditions()` within the `Fem1Module.cc`:

```c++
void Fem1Module::
_initBoundaryconditions()
{
  info() << "Init boundary conditions...";
    
  _applyOneBoundaryCondition("Cercle", 50.0);
  _applyOneBoundaryCondition("Bas", 5.0);
  _applyOneBoundaryCondition("Haut", 21.0);
}
```

the `_applyOneBoundaryCondition()` function is used to input the three Dirichlet boundary conditions here. For example, via  `_applyOneBoundaryCondition("Bas", 5.0);` we impose $5.0$ $^\text{o}\text{C}$ on the border "Bas" (french for bottom). Note, "Bas" is a physical tag in the `plancher.msh` file.  So in the snippet above, three Dirichlet conditions are applied ($50 ^\text{o}\text{C}, 5.0 ^\text{o}\text{C}, 21.0 ^\text{o}\text{C}$)  on three borders ('cercle', 'Bas', 'Haut').

### Post Process ###

For post processing the `ensight.case` file is outputted, which can be read by PARAVIS. The output is of the $\mathcal{P}_1$ FE order (on nodes).  
