
This document is intended to guide users in compiling ArcaneFEM on various clusters with different architectures. Installing codes on clusters can be challenging, so this serves as a starting point to simplify the process. The focus is on providing practical help rather than being overly detailed or perfect.

<details>
  <summary>
  
  # Adastra 
</summary>
In 2024, Adastra is the Fastest supercomputer in France. 
Adastra is a French supercomputer hosted at CINES, a Tier 1 computing site located in Montpellier. 
The Adastra supercomputer is an HPE-Cray EX system, combined with two ClusterStor E1000 storage systems.
Please have a look at this [link](https://dci.dci-gitlab.cines.fr/webextranet/architecture/index.html) to know more about Adastra architecture.
Here we shall concern ourselves with MI250X partition of Adastra, which are accelerated nodes of Adastra specialized for General Purpose computation 
on GPUs (GPGPU) computations (1 AMD Trento EPYC 7A53 64 cores 2.0 GHz processor with 256 Gio of DDR4-3200 MHz CPU memory per node, 4 Slingshot 200 Gb/s NICs, 
8 GPUs devices (4 AMD MI250X accelerator, each with 2 GPUs) with a total of 512 Gio of HBM2 per node).

In order to be able to run ArcaneFEM on Adastra in multi-gpu mode, you will need to compile dotnet, Arcane, and then ArcaneFEM, in that order. An example compilation is explained below.
  

  <details>
    <summary>

  #### Compile dotnet on Adastra
  </summary>
  
  Arcane depends on `.Net`, this dependencie is not found on Adastra, hence we begin by installaing dotnet. 
  If you have this already installed please skip this section.

  - On your personal PC, Locally download dotNet `tar.gz` and move it to move it to Adastra.
    Note please fill in your `ADASTRA_USERNAME`, and your `ADASTRA_PROJECT`

```bash
export ADASTRA_USERNAME="XXXXX"
export ADASTRA_PROJECT="YYYYY"

wget https://download.visualstudio.microsoft.com/download/pr/db901b0a-3144-4d07-b8ab-6e7a43e7a791/4d9d1b39b879ad969c6c0ceb6d052381/dotnet-sdk-8.0.401-linux-x64.tar.gz

scp dotnet-sdk-8.0.401-linux-x64.tar.gz ${ADASTRA_USERNAME}@adastra2.cines.fr:/lus/work/RES1/${ADASTRA_PROJECT}/${ADASTRA_USERNAME}/.
```

  - On Adastra create working directory

```bash
export ROOT_DIR=${WORKDIR}/ArcaneFEM
mkdir -p $ROOT_DIR
mkdir -p $ROOT_DIR/dotnet
cd $ROOT_DIR
```
- Extract the archive on Adastra

```bash
cd $ROOT_DIR/dotnet
tar xvzf ./../../dotnet-sdk-8.0.401-linux-x64.tar.gz
rm -rf ./../../dotnet-sdk-8.0.401-linux-x64.tar.gz
```
 </details>
 
  <details>
    <summary>

  #### Compile Arcane with GPU support on Adastra
  </summary>
 
 - Make sure `dotnet` is installed as it is a prerequiste of Arcane, see section above on dotnet installation. If it already installed you can proceed.
 - In the following you will need to add your assocatied allocation project in the `ADASTRA_PROJECT` variable.

  ```bash
export ADASTRA_PROJECT="XXXX"

export ROOT_DIR=${WORKDIR}/ArcaneFEM
mkdir -p $ROOT_DIR
cd $ROOT_DIR

cd $ROOT_DIR && git clone https://github.com/arcaneframework/framework.git && cd framework && git submodule update --init --recursive

salloc -A ${ADASTRA_PROJECT} --constraint=MI250 --nodes=1 -c 32  --gpus-per-node=1 --time=1:00:00

module purge

module load craype-x86-trento craype-accel-amd-gfx90a PrgEnv-gnu/8.5.0 rocm/5.7.1 GCC-GPU-3.0.0 hypre/2.29.0-mpi-omp  boost/1.83.0-mpi-python cray-hdf5-parallel/1.12.2.9 cmake/3.27.7

export ROOT_DIR=${WORKDIR}/ArcaneFEM
export SOURCE_DIR=${ROOT_DIR}/framework
export BUILD_DIR=${ROOT_DIR}/framework-build-release
export INSTALL_DIR=${ROOT_DIR}/framework-install

export DOTNET_ROOT=${ROOT_DIR}/dotnet
export PATH=${ROOT_DIR}/dotnet:${PATH}

export PREFIX_PATH="${ROCM_PATH};${ROCM_PATH}/hip;/opt/software/USERS_SOFTWARES/TRIOU/TRUST/TRUST-1.9.4-gpu/lib/src/LIBPETSC/petsc/linux_opt"
export HIP_ARCHITECTURES=gfx90a

export CC="${ROCM_PATH}/bin/amdclang"
export CXX="${ROCM_PATH}/bin/amdclang++"

export HDF5_CC="amdclang"
export HDF5_CLINKER="amdclang"
export HDF5_CXX="amdclang++"
export HDF5_CXXLINKER="amdclang++"

cmake -S ${SOURCE_DIR} -B ${BUILD_DIR} -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_PREFIX_PATH=${PREFIX_PATH} -DARCANE_ACCELERATOR_MODE=ROCMHIP -DCMAKE_HIP_ARCHITECTURES=${HIP_ARCHITECTURES} -DARCCORE_CXX_STANDARD=20 -DCMAKE_BUILD_TYPE=Release -DCMAKE_DISABLE_FIND_PACKAGE_SWIG=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_PETSc=TRUE -D_HYPRE_CONFIG_PATH=/opt/software/gaia/prod/3.0.0/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_plac/hypre-2.29.0-gcc-12.1-oah4/include/HYPRE_config.h

cd $BUILD_DIR
make -j all install
  ```
**NOTE:** you will need to patch a file `ArcaneMpi.cc`    replace `is_aware = (MPIX_Query_hip_support()==1);`  to `is_aware = 0;`
</details>

  <details>
    <summary>
  
  #### Compile ArcaneFEM with GPU support on Adastra
  </summary>
 - In the following you will need to add your assocatied allocation project in the `ADASTRA_PROJECT` variable.

  ```bash
export ADASTRA_PROJECT="XXXX"

export ROOT_DIR=${WORKDIR}/ArcaneFEM
mkdir -p $ROOT_DIR
cd $ROOT_DIR

cd $ROOT_DIR && git clone https://github.com/arcaneframework/arcanefem.git

salloc -A ${ADASTRA_PROJECT} --constraint=MI250 --nodes=1 -c 32  --gpus-per-node=1 --time=1:00:00

module purge

module load craype-x86-trento craype-accel-amd-gfx90a PrgEnv-gnu/8.5.0 rocm/5.7.1 GCC-GPU-3.0.0 hypre/2.29.0-mpi-omp  boost/1.83.0-mpi-python cray-hdf5-parallel/1.12.2.9 cmake/3.27.7


export PREFIX_PATH="${ROCM_PATH};${ROCM_PATH}/hip;/opt/software/USERS_SOFTWARES/TRIOU/TRUST/TRUST-1.9.4-gpu/lib/src/LIBPETSC/petsc/linux_opt"
export HIP_ARCHITECTURES=gfx90a

export CC="${ROCM_PATH}/bin/amdclang"
export CXX="${ROCM_PATH}/bin/amdclang++"

export HDF5_CC="amdclang"
export HDF5_CLINKER="amdclang"
export HDF5_CXX="amdclang++"
export HDF5_CXXLINKER="amdclang++"

export ROOT_DIR=${WORKDIR}/ArcaneFEM
export ARCANE_INSTALL_DIR=${ROOT_DIR}/framework-install
export BUILD_DIR=${ROOT_DIR}/arcanefem-build-release
export SOURCE_PATH=${ROOT_DIR}/arcanefem

cmake -S ${SOURCE_PATH} -B ${BUILD_DIR} -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_DIR} -DCMAKE_EXE_LINKER_FLAGS="-pthread" -DCMAKE_HIP_ARCHITECTURES=${HIP_ARCHITECTURES} -DHypre_INCLUDE_DIRS=/opt/software/gaia/prod/3.0.0/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_plac/hypre-2.29.0-gcc-12.1-oah4/include -DHypre_LIBRARY=/opt/software/gaia/prod/3.0.0/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_placeholder__/__spack_path_plac/hypre-2.29.0-gcc-12.1-oah4/lib

cd $BUILD_DIR
make -j all install
  ```
</details>
</details>

<details>
  <summary>
  
  # Topaz 
</summary>
Topaz is a French supercomputer hosted at CEA, designed for both CPU and GPU computing. It features 994 CPU nodes, each with 2 AMD EPYC Milan 7763 processors (128 cores per node) and 256 GB of memory, delivering 4.5 Pflops peak performance with InfiniBand HDR-100 interconnects. Additionally, it includes 75 hybrid nodes, each with 2 AMD EPYC Milan 7763 processors, 4 NVIDIA A100 GPUs, 512 GB of memory, and InfiniBand HDR interconnects, achieving 4.3 Pflops peak performance.

To run ArcaneFEM on Topaz, you need to compile dotnet, Arcane, and ArcaneFEM, following the steps detailed below.

<details>
    <summary>

  #### Compile dotnet on Topaz
  </summary>
  
Arcane depends on `.Net`, this dependencie is not found on Topaz, hence we begin by installaing dotnet. If you have this already installed please skip this section.

  - On your personal PC, locally download dotNet `tar.gz` and move it to move it to Topaz. Note please fill in your `TOPAZ_USERNAME`, and your `TOPAZ_PROJECT`

```bash
export TOPAZ_USERNAME="XXXXX"
export TOPAZ_PROJECT="contXXX/XXX"

wget https://download.visualstudio.microsoft.com/download/pr/db901b0a-3144-4d07-b8ab-6e7a43e7a791/4d9d1b39b879ad969c6c0ceb6d052381/dotnet-sdk-8.0.401-linux-x64.tar.gz

scp dotnet-sdk-8.0.401-linux-x64.tar.gz ${TOPAZ_USERNAME}@topaze.ccc.cea.fr:/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/.
```

  - On Topaz create working directory

```bash
export ROOT_DIR=${CCCWORKDIR}/ArcaneFEM
mkdir -p $ROOT_DIR
mkdir -p $ROOT_DIR/dotnet
cd $ROOT_DIR
```
- Extract the archive on Topaz

```bash
cd $ROOT_DIR/dotnet
tar xvzf ./../../dotnet-sdk-8.0.401-linux-x64.tar.gz
rm -rf ./../../dotnet-sdk-8.0.401-linux-x64.tar.gz
```
</details>

<details>
    <summary>

  #### Compile Arcane on Topaz
  </summary>
  
  - On your personal PC, locally checkout Arcane from GitHub and move it to Topaz. Note please fill in your `TOPAZ_USERNAME`, and your `TOPAZ_PROJECT`

```bash
export TOPAZ_USERNAME="XXXXX"
export TOPAZ_PROJECT="contXXX/XXX"

git clone https://github.com/arcaneframework/framework.git && cd framework && git submodule update --init --recursive && cd ..

scp -r framework ${TOPAZ_USERNAME}@topaze.ccc.cea.fr:/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/.

rm -rf framework
```

- On Topaz you will need to do the following. Please note, here the idea is to pick up `hypre` and `parmetis` that was precopiled and add it locally on Topaz, since the locally avaliable packages on Topaz are not compile with CUDA version 12 support on Topaz and this is necessary for C++ 20 that ArcaneFEM absolutely needs.
  - move framework to `ArcaneFEM` directory

```bash
export ROOT_DIR=${CCCWORKDIR}/ArcaneFEM
mkdir -p $ROOT_DIR
cd $ROOT_DIR
mv ${CCCWORKDIR}/framework .
```
  - configure Arcane framework

```bash
export TOPAZ_USERNAME="XXXXX"
export TOPAZ_PROJECT="contXXX/XXX"

export ROOT_DIR=${CCCWORKDIR}/ArcaneFEM
mkdir -p $ROOT_DIR
cd $ROOT_DIR

module purge

module load gnu/11 flavor/openmpi/cuda-12.2 mpi/openmpi/4.1.5 cmake/3.26 hdf5/1.12.0 metis/5.1.0

export ROOT_DIR=${CCCWORKDIR}/ArcaneFEM
export SOURCE_DIR=${ROOT_DIR}/framework
export BUILD_DIR=${ROOT_DIR}/framework-build-release
export INSTALL_DIR=${ROOT_DIR}/framework-install

export PATH=/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/ArcaneFEM/dotnet:$PATH

_HWLOC_PATH="/ccc/products/hwloc-2.2.0/system/default"
_OTF2_PATH="/ccc/products/otf2-2.3/nvidia--22.2__openmpi--4.0.1/default"
_PETSC_PATH="/ccc/products/petsc-3.17.4/gcc--11.1.0__openmpi--4.0.1/default"

HYPRE_PATH="/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/ArcaneFEM/hypre/2.29.0/gpu"
PARMETIS_PATH="/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/ArcaneFEM/parmetis/4.0.3-ompi405"

COMMON_CMAKE_PREFIX_PATH="${_HWLOC_PATH};${_OTF2_PATH};${_PETSC_PATH};${HYPRE_PATH};${PARMETIS_PATH}"

export CXX=`which c++`
export CC=`which gcc`
export CXX CC

cmake --preset Arcane -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane -DCMAKE_PREFIX_PATH="${COMMON_CMAKE_PREFIX_PATH}" -DARCANE_WANT_CUDA=TRUE -S ${SOURCE_DIR} -B ${BUILD_DIR} -DARCANE_CUSTOM_MPI_DRIVER="/usr/bin/ccc_mprun" -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_BUILD_TYPE=Release -DARCCORE_CXX_STANDARD=20  -DCMAKE_DISABLE_FIND_PACKAGE_SWIG=TRUE
```

- finally make and make install

```bash
cd $BUILD_DIR
make -j all install
```
  </details>

<details>
    <summary>

  #### Compile ArcaneFEM on Topaz
  </summary>

  - On your personal PC, locally checkout ArcaneFEM from GitHub and move it to Topaz. Note please fill in your `TOPAZ_USERNAME`, and your `TOPAZ_PROJECT`

```bash
export TOPAZ_USERNAME="XXXXX"
export TOPAZ_PROJECT="contXXX/XXX"

git clone https://github.com/arcaneframework/arcanefem.git

scp -r arcanefem ${TOPAZ_USERNAME}@topaze.ccc.cea.fr:/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/.

rm -rf arcanefem
```

- On Topaz you will need to do the following.
    - move arcanefem to `ArcaneFEM` directory

```bash
export ROOT_DIR=${CCCWORKDIR}/ArcaneFEM
mkdir -p $ROOT_DIR
cd $ROOT_DIR
mv ${CCCWORKDIR}/arcanefem .
```
  - configure ArcaneFEM framework

```bash
export TOPAZ_USERNAME="XXXXX"
export TOPAZ_PROJECT="contXXX/XXX"

export ROOT_DIR=${CCCWORKDIR}/ArcaneFEM
mkdir -p $ROOT_DIR
cd $ROOT_DIR

module purge
module load gnu/11 flavor/openmpi/cuda-12.2 mpi/openmpi/4.1.5 cmake/3.26 hdf5/1.12.0 metis/5.1.0

export ROOT_DIR=${CCCWORKDIR}/ArcaneFEM
export SOURCE_DIR=${ROOT_DIR}/arcanefem
export BUILD_DIR=${ROOT_DIR}/arcanefem-build-release
export ARCANE_DIR=${ROOT_DIR}/framework-install

export PATH=/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/ArcaneFEM/dotnet:$PATH

_HWLOC_PATH="/ccc/products/hwloc-2.2.0/system/default"
_OTF2_PATH="/ccc/products/otf2-2.3/nvidia--22.2__openmpi--4.0.1/default"
_PETSC_PATH="/ccc/products/petsc-3.17.4/gcc--11.1.0__openmpi--4.0.1/default"

HYPRE_PATH="/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/ArcaneFEM/hypre/2.29.0/gpu"
PARMETIS_PATH="/ccc/work/${TOPAZ_PROJECT}/${TOPAZ_USERNAME}/ArcaneFEM/parmetis/4.0.3-ompi405"

COMMON_CMAKE_PREFIX_PATH="${_HWLOC_PATH};${_OTF2_PATH};${_PETSC_PATH};${HYPRE_PATH};${PARMETIS_PATH};${ARCANE_DIR}"

cmake -S ${SOURCE_DIR} -B ${BUILD_DIR}  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=${COMMON_CMAKE_PREFIX_PATH}
```

- finally make and make install

```bash
cd $BUILD_DIR
make -j all install
```

</details>

</details>

