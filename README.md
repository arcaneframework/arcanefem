Very simple codes to test Finite Element Methods using Arcane

You need to install a recent (3.7+) version of Arcane Framework before using this sample.

#### How to run the solver  ####

To compile the sources (in a directory)

~~~{sh}
ARCANE_INSTALL_DIR=/path/to/arcane/installation
BUILD_DIR=/tmp/build
SOURCE_PATH=/path/to/sources
cmake -S ${SOURCE_PATH} -B ${BUILD_DIR} -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_DIR}
cmake --build ${BUILD_DIR}
~~~

To execute an example from poisson solver (in a directory)

~~~{sh}
cd ${BUILD_DIR}/poisson && ./FemTest Test.conduction.arc
~~~

After running the test case, you can display the results with ParaView:

~~~bash
paraview ${BUILD_DIR}/poisson/output/depouillement/ensight.case
~~~

## Todo List ##

#### Short term ####
- [ ] Row-Column elimination for Dirichlet BC
- [ ] Convection operator for transient heat conduction problem
- [x] Transient Heat conduction problem
- [x] Update documentation
- [x] Point-Dirichlet for elasticity
- [x] Traction condition for elasticity
- [x] Possibility to change linear-algebra backed in .arc files
- [x] Reading group of nodes (Arcane::NodeGroup) from msh file 
- [x] Point Dirichlet boundary condition
- [x] First vectorial FEM example - Bilaplacian
- [x] Reorganization (FemLib, Poission, Bilaplacian, .... )

#### Long term ####

- [ ] Parallel visualization
- [ ] top-ii-vol integration
- [ ] FEM interpolation
- [ ] MED file testing
