# Installation

## Prerequisites

Before installing ArcaneFEM, ensure you have the following dependencies installed on your system:

#### Required Dependencies

- **Arcane Framework** (version ≥ 3.14.14)
- **CMake** (version ≥ 3.18)
- **C++ Compiler** with C++20 support (GCC ≥ 10, Clang ≥ 12, or equivalent)
- **MPI** implementation (OpenMPI, MPICH, or Intel MPI)

#### Optional Dependencies

<!-- tabs:start -->

##### **NVIDIA / CUDA**

- **CUDA Toolkit** (for GPU acceleration)  
- **ParaView** (version > 5.12, for visualization)  
- **HDF5** (for advanced output formats)

##### **AMD ROCm**

- **ROCm Toolkit** (for GPU acceleration)  
- **ParaView** (version > 5.12, for visualization)  
- **HDF5** (for advanced output formats)

<!-- tabs:end -->


> **Note**: information on [how to build Arcane can be found here](https://arcaneframework.github.io/arcane/userdoc/html/d7/d94/arcanedoc_build_install.html)

## Installing ArcaneFEM

#### Step 1: Clone the Repository
```bash
git clone https://github.com/arcaneframework/arcanefem.git
cd arcanefem
```

#### Step 2: Set Environment Variables
```bash
# Set the path to your Arcane installation
export ARCANE_INSTALL_DIR=/path/to/arcane/installation

# Optional: Add to your shell profile (.bashrc, .zshrc, etc.)
echo "export ARCANE_INSTALL_DIR=/path/to/arcane/installation" >> ~/.bashrc
```
> Note: How to compile Arcane has been explained in the next section. 

#### Step 3: Configure Build
```bash
# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. -DCMAKE_PREFIX_PATH=${ARCANE_INSTALL_DIR} \
         -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_CXX_STANDARD=20
```

#### Step 4: Build
```bash
# Build the project
cmake --build . --parallel $(nproc)

# Or using make directly
make -j$(nproc)
```

#### Step 5: Run tests
```bash
# Build the project
ctest
```

#### Step 6: Verify Installation
```bash
# Test the elasticity solver
cd elasticity

# Run a simple test case
./Elasticity Test.Elasticity.arc
```

## System-Specific Install

<!-- tabs:start -->

### **Ubuntu 24.04 NVIDIA**
- Install required packages

```bash
sudo apt-get update
sudo apt install -y apt-utils build-essential iputils-ping python3 swig \
  git g++-12 gcc-12 gfortran nvidia-cuda-toolkit \
  libglib2.0-dev libxml2-dev libhdf5-openmpi-dev \
  libparmetis-dev libunwind-dev dotnet8 cmake \
  libhypre-dev libpetsc-real-dev libtrilinos-teuchos-dev libtrilinos-epetra-dev \
  libtrilinos-tpetra-dev libtrilinos-kokkos-dev libtrilinos-ifpack2-dev \
  libtrilinos-ifpack-dev libtrilinos-amesos-dev libtrilinos-galeri-dev \
  libtrilinos-xpetra-dev libtrilinos-epetraext-dev \
  libtrilinos-triutils-dev libtrilinos-thyra-dev \
  libtrilinos-kokkos-kernels-dev libtrilinos-rtop-dev \
  libtrilinos-isorropia-dev libtrilinos-belos-dev
```

- Cmake configure Arcane Framework

```bash
cmake -S /path/to/Arcane/sources -B /path/to/Arcane/build \
   -DCMAKE_INSTALL_PREFIX=$HOME/framework-install \
   -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_COMPILER=g++-12 \
   -DCMAKE_C_COMPILER=gcc-12 \
   -DARCANE_ACCELERATOR_MODE=CUDANVCC \
   -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc \
   -DARCCORE_CXX_STANDARD=20
```

- Install Arcane Framework

```bash
cd /path/to/Arcane/build
make -j all
make install
```

- Cmake configure ArcaneFEM

```bash
cmake -S /path/to/ArcaneFEM/sources -B /path/to/ArcaneFEM/build \ 
   -DCMAKE_PREFIX_PATH=$HOME/framework-install \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_STANDARD=20
```

- Install ArcaneFEM

```bash
cd /path/to/ArcaneFEM/build
make -j all
```

- Testing

```bash
ctest
```

### **Ubuntu 24.04 AMD**
- Install required packages

```bash
sudo apt-get update
sudo apt install -y apt-utils build-essential iputils-ping python3 swig \
  git g++-12 gcc-12 gfortran rocm amdgpu-dkms \
  libglib2.0-dev libxml2-dev libhdf5-openmpi-dev \
  libparmetis-dev libunwind-dev dotnet8 cmake \
  libhypre-dev libpetsc-real-dev libtrilinos-teuchos-dev libtrilinos-epetra-dev \
  libtrilinos-tpetra-dev libtrilinos-kokkos-dev libtrilinos-ifpack2-dev \
  libtrilinos-ifpack-dev libtrilinos-amesos-dev libtrilinos-galeri-dev \
  libtrilinos-xpetra-dev libtrilinos-epetraext-dev \
  libtrilinos-triutils-dev libtrilinos-thyra-dev \
  libtrilinos-kokkos-kernels-dev libtrilinos-rtop-dev \
  libtrilinos-isorropia-dev libtrilinos-belos-dev
```
> Note: Follow ROCm and AMDGPU [installation instructions](https://rocm.docs.amd.com/projects/install-on-linux/en/latest/install/quick-start.html) to have the `rocm` and `amdgpu-dkms` packages avaliable for install.


- Set env variables

```bash
export ROCM_ROOT=/opt/rocm-6.4.1-60401
export CC=/opt/rocm/llvm/bin/clang
export CXX=/opt/rocm/llvm/bin/clang++
export CMAKE_HIP_COMPILER=/opt/rocm/hip/bin/hipcc
```

- Cmake configure Arcane Framework (use rocminfo to get DCMAKE_HIP_ARCHITECTURES)

```bash
cmake -S /path/to/Arcane/sources -B /path/to/Arcane/build \
   -DCMAKE_INSTALL_PREFIX=$HOME/framework-install \
   -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane \
   -DCMAKE_BUILD_TYPE=Release \
   -DARCANE_ACCELERATOR_MODE=ROCMHIP \
   -DCMAKE_HIP_ARCHITECTURES=gfx90a \
   -DARCCORE_CXX_STANDARD=20
```

- Install Arcane Framework

```bash
cd /path/to/Arcane/build
make -j all
make install
```

- Cmake configure ArcaneFEM

```bash
cmake -S /path/to/ArcaneFEM/sources -B /path/to/ArcaneFEM/build \ 
   -DCMAKE_PREFIX_PATH=$HOME/framework-install \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_STANDARD=20
```

- Install ArcaneFEM

```bash
cd /path/to/ArcaneFEM/build
make -j all
```

- Testing

```bash
ctest
```

### **macOS 15.5**

- Using Homebrew

```bash
brew install git cmake gcc open-mpi hdf5-mpi scotch glib dotnet
```

- Cmake configure Arcane Framework

```bash
cmake -S /path/to/Arcane/sources -B /path/to/Arcane/build \
    -DCMAKE_INSTALL_PREFIX=$HOME/framework-install \
    -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane \
    -DCMAKE_BUILD_TYPE=Release \
    -DARCCORE_CXX_STANDARD=20 \
    -DCMAKE_CXX_COMPILER=g++-15 \
    -DCMAKE_C_COMPILER=gcc-15 \
    -DCMAKE_DISABLE_FIND_PACKAGE_SWIG=TRUE \
    -DCMAKE_DISABLE_FIND_PACKAGE_Hypre=TRUE \
    -DCMAKE_DISABLE_FIND_PACKAGE_PETSc=TRUE
```

- Install Arcane Framework

```bash
cd /path/to/Arcane/build
make -j$(sysctl -n hw.ncpu)
make install
```

- Cmake configure ArcaneFEM

```bash
cmake -S /path/to/ArcaneFEM/sources -B /path/to/ArcaneFEM/build \ 
   -DCMAKE_PREFIX_PATH=$HOME/framework-install \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_STANDARD=20
```

- Install ArcaneFEM

```bash
cd /path/to/ArcaneFEM/build
make -j$(sysctl -n hw.ncpu)
```

- Testing

```bash
ctest
```

<!-- tabs:end -->





## Troubleshooting

### Common Issues

1. **Arcane Framework Not Found**
   ```bash
   # Ensure ARCANE_INSTALL_DIR is correctly set
   export ARCANE_INSTALL_DIR=/correct/path/to/arcane
   ```

2. **MPI Compilation Errors**
   ```bash
   # Load MPI module or set compiler wrappers
   export CC=mpicc
   export CXX=mpicxx
   ```

3. **CUDA Compilation Issues**
   ```bash
   # Ensure CUDA toolkit is in PATH
   export PATH=/usr/local/cuda/bin:$PATH
   export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
   ```

### Getting Help

- Check the [Arcane Framework Documentation](https://arcaneframework.github.io/arcane/userdoc/html/index.html)
- Open an issue on the [ArcaneFEM GitHub repository](https://github.com/arcaneframework/arcanefem/issues)