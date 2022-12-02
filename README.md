Very simple code to test Finite Element Methods using Arcane

You need to install a recent (3.7+) version of Arcane Framework before using this sample.
To compile and execute (in direcrot

~~~{sh}
ARCANE_INSTALL_DIR=/path/to/arcane/installation
BUILD_DIR=/tmp/build
SOURCE_PATH=/path/to/sources
cmake -S ${SOURCE_PATH}/fem1 -B ${BUILD_DIR} -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_DIR}
cmake --build ${BUILD_DIR}
cd ${BUILD_DIR} && ./FemTest1 FemTest1.arc
~~~

After running the test case, you can display the results with ParaView:

~~~
paraview ${BUILD_DIR}/output/depouillement/ensight.case
~~~
