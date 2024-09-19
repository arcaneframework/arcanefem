Welcome to the repository showcasing Finite Element Method (FEM) based solvers developed using the Arcane Framework. The FEM solvers/algorithms here are optimized for both CPU and GPU-based parallel computing environments.

Before diving into the samples provided, please ensure you have installed a recent version (3.14.4) of the Arcane Framework.

## How to test  ##

It is simple **compile**$\rightarrow$**execute**$\rightarrow$**visualize**


### compile ###
To compile the sources, follow these steps:

~~~{sh}
# Set up paths
ARCANE_INSTALL_DIR=/path/to/arcane/installation
BUILD_DIR=/tmp/build
SOURCE_PATH=/path/to/sources

# Invoke CMake to configure the build
cmake -S ${SOURCE_PATH} -B ${BUILD_DIR} -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_DIR}

# Build the project
cmake --build ${BUILD_DIR}
~~~

### Execute ###

Once compiled, execute any module of your choice. For example for the elasticity solver. Navigate to the appropriate directory: 

~~~{sh}
cd ${BUILD_DIR}/elasticity
~~~
Then run the executable with the desired input file:
~~~{sh}
./Elasticity Test.Elasticity.arc
~~~

for parallel run (domain-decompostion) e.g:

```sh
mpirun -n 4 ./Elasticity Test.Elasticity.arc
```

Alternatively, you can provide command-line arguments to run the solver:

~~~{sh}
./Elasticity -A,CaseDatasetFileName=Test.Elasticity.arc
~~~
For additional commands to control Arcane, refer to  [Arcane Documentation](https://arcaneframework.github.io/arcane/userdoc/html/d8/dd6/arcanedoc_execution_launcher.html).* 

### Visualize ###

After running the test cases visualize the results using ParaView:

~~~{sh}
paraview ${BUILD_DIR}/elastcity/output/depouillement/vtkhdfv2/Mesh0.hdf
~~~

please note you will need the latest ParaView ( > 5.12) and Arcane framework compiled with mpi support for hdf5. 
