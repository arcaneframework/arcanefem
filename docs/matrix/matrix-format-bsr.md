# Block Sparse Row Format Documentation


`BSRFormat` class provides facilities for storing, manipulating and assembling the global stiffness matrix from mesh informations.

- It supports both triangular and tetrahedral elements with multiple degrees-of-freedom.
  
- It supports `Aleph` and `Hypre` linear system backend.
  
- It is GPU-accelerated (enabled with `-A,AcceleratorRuntime` option).
  
- It uses a `BSRMatrix` to represent the global stiffness matrix under the hood. `BSRMatrix` stores the non-zero coefficients of a matrix within three arrays: `values`, `columns` and `row_index` as explained [here](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2024-0/sparse-blas-bsr-matrix-storage-format.html).

ArcaneFEM's `Elasticity` and `Testlab` modules rely on it.

## Initializing BSRFormat in your Module

```cpp
#include "BSRFormat.h"
#include <arcane/accelerator/core/IAcceleratorMng.h>
#include "arcane/accelerator/VariableViews.h"
```

To use `BSRFormat` in your module, it is recommended to include it as a member field of the module.

Ensure you initialize it by calling its constructor:

```cpp
BSRFormat(ITraceMng* tm, RunQueue& queue, const FemDoFsOnNodes& dofs_on_nodes)
```

Example:

```cpp
class FemModule : public ArcaneFemObject {

  FemModule(const ModuleBuildInfo& mbi)
  : ArcaneFemObject(mbi),
  // ...
  , m_dofs_on_nodes(mbi.subDomain()->traceMng())
  , m_bsr_format(mbi.subDomain()->traceMng(), *(mbi.subDomain()->acceleratorMng()->defaultQueue()), m_dofs_on_nodes)    // Add this line for BSR
  { ... }
  // ...
  FemDoFsOnNodes m_dofs_on_nodes;
  BSRFormat m_bsr_format;                // Add this line for BSR
}
```

After constructing the `BSRFormat` object, you must initialize it using the `initialize(...)` method:

```cpp
// Prototype:
void initialize(IMesh* mesh, Int32 nb_edge, bool does_linear_system_use_csr);

// Parameters:
//   - mesh: The mesh object to associate with the `BSRFormat`.
//   - nb_edge: In 3D, this represents the number of edges; in 2D, it represents the number of faces.
//   - does_linear_system_use_csr: Set to `true` if Hypre is used as the backend (CSR format), otherwise `false`.
```

Implementation example of the initialization process in the `startInit` method:

```cpp
void FemModule::startInit() {
  // ...
  if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
    auto use_csr_in_linear_system = options()->linearSystem.serviceName() == "HypreLinearSystem";
    if (m_matrix_format == "BSR")
      m_bsr_format.initialize(mesh(), 1, use_csr_in_linear_system, 0);
    else
      m_bsr_format.initialize(mesh(), 1, use_csr_in_linear_system, 1);
    m_bsr_format.computeSparsity();
  }
}
```

Note that here are are adding the two types of BSR available in ArcaneFEM, the atomic-free one and the normal BSR. If you just want one of the two adapt the code above. Also note that here we computed the sparsity of the matrix `m_bsr_format.computeSparsity();`  as we imagine the sparsity is constant. 

## Bilinear Operator Assembly

`BSRFormat` offers two approaches for assembling the bilinear operator i.e. global stiffness matrix:
1. Default Approach (`m_matrix_format == "BSR"`)
2. Atomic-free Approach (`m_matrix_format == "AF-BSR"`)

Depending on the method chosen, the assembly process now has the task:

- either to assemble the **local element matrix** contributions and  accumulate these contributions of into the global matrix. This will be the case for *Default Approach*
- or to assemble the **local element vector** contributions and  accumulate these contributions of into the global matrix. This will be the case for *Atomic-free Approach*



<details open>
  <summary><h3>Default Approach</h3></summary>

The default way of assembling the global stiffness matrix is by integrating over the elements of the mesh.
Because of that and of GPU compatibility (i.e. shared-memory parallelism), this method uses an `atomic` operation within its implementation.

> In theory, these operations can lead to non-deterministic behavior, i.e. to a matrix whose coefficients can change from one run to the next for the same simulation.
> In practice, no such behavior has been observed by us, even in large test case.

#### Contributions

```cpp
#include <arcane/accelerator/VariableViews.h>

namespace ax = Arcane::Accelerator;

// ...

ARCCORE_HOST_DEVICE RealMatrix<3, 3> computeElementMatrixTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord) { ... }

// ...

ARCCORE_HOST_DEVICE RealMatrix<1, 3> computeElementVectorTria3Gpu(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord, Int32 node_lid)
{ ... }

// ...

void FemModule::_assembleBilinearOperator() {
  // ...
    UnstructuredMeshConnectivityView m_connectivity_view(mesh());
    auto cn_cv = m_connectivity_view.cellNode();
    auto m_queue = subDomain()->acceleratorMng()->defaultQueue();
    auto command = makeCommand(m_queue);
    auto in_node_coord = ax::viewIn(command, m_node_coord);
    m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord); });

  // ...
}
```

</details>

<details open>
  <summary><h3>Atomic-free Approach</h3></summary>

This method assembles the global stiffness by integrating over the DOF's  (nodes since we are doing   $P_1$  FEM) of the mesh. It does not use `atomic` operations and therefore is deterministic.

#### Contributions Assembly

```cpp
#include <arcane/accelerator/VariableViews.h>

namespace ax = Arcane::Accelerator;

// ...

ARCCORE_HOST_DEVICE RealMatrix<3, 3> computeElementMatrixTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord) { ... }

void FemModule::_assembleBilinearOperator() {
  // ...
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto queue = subDomain()->acceleratorMng()->defaultQueue();
  auto command = makeCommand(queue);
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, node_lid); });
}
```

</details>

*NOTE* Additional thing to do is to make sure that the lambda functions do not have private or protected access within its class. In the example above `_assembleBilinearOperator` needs to be `public`, this is only if you compile using `nvcc` if `clang` is used this step is not needed.

### Translate BSRFormat to DoFLinearSystem

ArcaneFEM uses `DoFLinearSystem` class to represent and handle linear systems. You can translate the BSRFormat to a given linear system using the `toLinearSystem(DoFLinearSystem &ls)` method.

Normally at this point your matrix assembly phase is on GPUs. You could for instance use matrix-assembly on GPU and matrix solving on CPU's. A `m_bsr_format.toLinearSystem(m_linear_system);` is required for that.  You would additionally need to initialize the linear system in a proper way
```c++
m_linear_system.initialize(subDomain(), acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");
```



## Linear Operator Assembly 

The `setValue(DoFLocalId row, DoFLocalId col, Real value)` method of `BSRMatrix` can be used to set coefficients in the matrix.
You can access the representation of `BSRFormat` i.e. `BSRMatrix` by using the `matrix()` method.

Example (applying Dirichlet boundary condition):
```cpp
// ...
Real Penalty = options()->penalty(); // 1.0e30 is the default
ENUMERATE_ (Node, inode, ownNodes()) {
  NodeLocalId node_id = *inode;
  if (m_u_dirichlet[node_id]) {
    DoFLocalId dof_id = node_dof.dofId(*inode, 0);
    m_bsr_format.matrix().setValue(dof_id, dof_id, Penalty);
    Real u_g = Penalty * m_u[node_id];
    rhs_values[dof_id] = u_g;
  }
} 
```



<details>
  <summary><h2>Algorithms and Implementation Details</h2></summary>

  This part is dedicated to the implementation details of `BSRFormat`, in particular the bilinear assembly algorithms.

<details>
  <summary><h3>Default Bilinear Operator Assembly</h3></summary>

  #### Sparsity Construction Algorithm

  ##### 1. Populate `row_index` Array

  1. Compute the `neighbors` array: At index `i`, `neighbors` contains the number of neighbors i.e. the number of connected nodes of node `i`.

     a. Loop over the elements of the mesh in parallel. Store each edge of the element in the `edge` array. Edge `i` (`i: 0 -> nb_edge_per_element`) is stored at index `cur_element_idx * nb_edge_per_element + i`. Edges are represented using a 64-bit integer. The first 32 bits store the `id` of the smaller node in the edge, while the last 32 bits store the `id` of the larger node.

     b. Sort the `edges` array into `sorted_edges`.

     c. Loop over the `sorted_edges` array in parallel. For each edge, if `sorted_edges[cur_edge_idx + 1] != cur_edge`, increment `neighbors[src]` and `neighbors[dst]` by `1` with an atomic operation. This conditional is needed to ensure that we don't count the same edge multiple time (for edges that are shared between multiple elements of the mesh).

  2. `row_index` is the [exclusive scan](https://en.wikipedia.org/wiki/Prefix_sum#Inclusive_and_exclusive_scans) of `neighbors`.

 ##### 2. Populate `columns` Array

 1. Compute `sorted_edges` (see 1.b)
 2. Loop over the edges in parallel. For each edge, if `sorted_edges[cur_edge_idx + 1] != cur_edge`, "register" the link `src -- dst` in `columns`: Get the start position of the row `src` in the matrix at `row_index[src]`. Get the current offset in the `src` row at `offsets[src]`. Put `dst` at `columns[start + offset]`. Increment `offsets[src]` by `1`.
    
#### Contributions Assembly Algorithm

1. Loop over the elements in parallel.

      a. Compute the local stiffness matrix of the element.

      b. Loop over the nodes of the current element. Loop over the nodes of the current element. Add the contribution `local_element_matrix[node1, node2]` into `bsr_matrix.values[node1, node2]` using an atomic operation.

</details>
<details>
<summary><h3>Atomic-free Bilinear Operator Assembly</h3></summary>

#### Sparsity Construction Algorithm

The only step in the previous sparsity construction method where an atomic operation is used is at 1.c i.e. for computing the `neighbors` array.

To build the sparsity without atomic, we compute the `neighbors` array using [Arcane's node-node connectivity](https://github.com/arcaneframework/framework/pull/1614). This connectivity is computed by Arcane on-demand. It uses node-edge connectivity under the hood and is not accelerated i.e. it will not use GPU to compute the connectivity.

The rest of the algorithm doesn't change apart from the iterations over the edges which are done using cell-edge in 3D (and cell-face in 2D) connectivities of Arcane.

#### Contributions Assembly Algorithm

1. Loop over the nodes in parallel. Loop over the elements of the node.

      a. Compute the local stiffness matrix of the element.

      b. Loop over the nodes of the element. Add the contribution `local_element_matrix[node1, node2]` into `bsr_matrix[node1, node2]`. Atomic is not needed here.
</details>
</details>
