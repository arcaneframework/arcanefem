Very simple codes to test Finite Element Methods using Arcane

You need to install a recent (3.7+) version of Arcane Framework before using this sample.

## How to test a solver  ##

It is simple **compile**$\rightarrow$**execute**$\rightarrow$**visualize**

- **Compile** the sources (in a directory)

~~~{sh}
ARCANE_INSTALL_DIR=/path/to/arcane/installation
BUILD_DIR=/tmp/build
SOURCE_PATH=/path/to/sources
cmake -S ${SOURCE_PATH} -B ${BUILD_DIR} -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_DIR}
cmake --build ${BUILD_DIR}
~~~

- Now you can **execute** an example from  elasticity solver (in a directory) `elastcity`
~~~{sh}
cd ${BUILD_DIR}/elastcity
~~~
~~~{sh}
./Elasticity Test.Elasticity.arc
~~~

*Note that you can also run the solver via the following command* 
~~~{sh}
./Elasticity -A,CaseDatasetFileName=Test.Elasticity.arc
~~~
  *Note: There are other commands to control Arcane these can be found [here](https://arcaneframework.github.io/arcane/userdoc/html/d8/dd6/arcanedoc_execution_launcher.html).* 

- After running the test case, you can **visualize** the results with ParaView:

~~~{sh}
paraview ${BUILD_DIR}/elastcity/output/depouillement/ensight.case
~~~

## Todo List ##

#### Short term ####
- [ ] Elastodynamics
- [ ] PETSc command-line paramerts or set parameters via a string
- [ ] Multi-mesh support
- [ ] Adhoc Matrix and Vector assemblies
- [x] Row elimination for Dirichlet BC
- [x] Row-Column elimination for Dirichlet BC
- [x] New parameter to switch the type of method for imposing Dirichlet BC
- [x] Gradient operator
- [x] Convection operator for transient heat conduction problem
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
- [ ] FEM interpolation via MedCoupling
- [ ] MED file testing
