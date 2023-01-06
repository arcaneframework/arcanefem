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

To execute an example from poission solver (in a directory)

~~~{sh}
cd ${BUILD_DIR}/poission && ./FemTest Test.conduction.arc
~~~

After running the test case, you can display the results with ParaView:

~~~bash
paraview ${BUILD_DIR}/poission/output/depouillement/ensight.case
~~~

## Todo List ##

#### Short term ####

- [ ] Update documentation
- [ ] Point-Dirichlet for elasticity
- [ ] Traction condition for elasticity
- [ ] Possibility to change linear-algebra backed in .arc files
- [ ] Row-Column elimination for Dirichlet BC
- [x] Reading group of nodes (Arcane::NodeGroup) from msh file 
- [x] Point Dirichlet boundary condition
- [x] First vectorial FEM example - Bilaplacian
- [x] Reorganization (FemLib, Poission, Bilaplacian, .... )

#### Long term ####

- [ ] Parallel visualization
- [ ] top-ii-vol integration
- [ ] FEM interpolation
- [ ] MED file testing
