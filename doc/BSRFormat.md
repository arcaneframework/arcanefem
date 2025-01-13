# BSRFormat Documentation

--- 
**TODO:**
- [ ] Update "Linear Operation Assembly" with `ArcaneFemFunctions`
- [ ] "Implementation Details" part
---

`BSRFormat` class provides facilities for storing, manipulating and assembling the global stiffness matrix. It uses a `BSRMatrix` to store the matrix under the hood. It is GPU compatible and will use accelerators if `-A,AcceleratorRuntime` option is set.

> `BSRMatrix` represents a matrix using three arrays: `row_index`, `columns` and `values` as explained [here](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2024-0/sparse-blas-bsr-matrix-storage-format.html). Size `n` of a block `n X n` is the number of degrees-of-freedom (DoFs) for the simulation.

<details>
  <summary><h2>Integration Guide</h2></summary>

### Setting up

To use `BSRFormat` in your module you will need to make it a member field. You must call its constructor.

```cpp
// proto: BSRFormat(ITraceMng* tm, RunQueue& queue, const FemDoFsOnNodes& dofs_on_nodes)

class FemModule : public ArcaneFemObject {
  FemModule(const ModuleBuildInfo& mbi) {   // ctor
    // ...
    m_bsr_format = BSRFormat<1>(mbi.subDomain()->traceMng(), *(mbi.subDomain()->acceleratorMng()->defaultQueue()), m_dofs_on_nodes);
  }
  // ...
  BSRFormat<1> m_bsr_format;                // 1 degree-of-freedom
}
```

Then you must initialize it using the `initialize(...)` method:

```cpp
// proto: void initialize(IMesh* mesh, Int32 nb_edge, bool does_linear_system_use_csr)
//   - nb_edge: number of edge in 3D, number of face in 2D
//   - does_linear_system_use_csr: True if Hypre is the backend, False otherwise

void FemModule::startInit() {
  // ...
  auto use_csr_in_linear_system = options()->linearSystem.serviceName() == "HypreLinearSystem";
  auto nb_edge = mesh->dimension() == 2 ? nbFace() : nbEdge();
  m_bsr_format.initialize(mesh, nb_edge, use_csr_in_linear_system); // initialize BSRFormat
}
```
### Bilinear Operator Assembly

`BSRFormat` offers two approaches for assembling the bilinear operator i.e. global stiffness matrix:
1. Default
2. Atomic-free

For both approaches, ArcaneFEM divides the assembly in 2 sub-steps:
1. Construction of the sparsity of the structure
2. Assembly of the contributions of local element matrix

<details open>
  <summary><h4>Default</h4></summary>

The default way of assembling the global stiffness matrix is by integrating over the elements of the mesh.
Because of that and of GPU compatibility (i.e. shared-memory parallelism), this method uses an `atomic` operation within its implementation.

> In theory, these operations can lead to non-deterministic behavior, i.e. to a matrix whose coefficients can change from one run to the next for the same simulation.
> In practice, no such behavior has been observed by us, even in large test case.

##### Construction of the Sparsity

```cpp
void FemModule::compute() {
  // ...
  m_bsr_format.computeSparsity()
}
```
  
##### Assembly of the Contributions

```cpp
namespace ax = Arcane::Accelerator;

// ...

ARCCORE_HOST_DEVICE FixedMatrix<3, 3> computeElementMatrixTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord) { ... }

void FemModule::compute() {
  // ...
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto command = makeCommand(m_queue);
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_bsr_format.assembleBilinear([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3(cell_lid, cn_cv, in_node_coord); });
}
```

</details>

<details open>
  <summary><h4>Atomic-free</h4></summary>

This method assembles the global stiffness by integrating over the nodes of the mesh. It does not use `atomic` operations and therefore is deterministic.

##### Construction of the Sparsity

```cpp
void FemModule::compute() {
  // ...
  m_bsr_format.computeSparsityAtomicFree()
}
```

##### Assembly of the Contributions

```cpp
namespace ax = Arcane::Accelerator;

// ...

ARCCORE_HOST_DEVICE FixedMatrix<3, 3> computeElementMatrixTria3(CellLocalId cell_lid, const IndexedCellNodeConnectivityView& cn_cv, const ax::VariableNodeReal3InView& in_node_coord) { ... }

void FemModule::compute() {
  // ...
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto command = makeCommand(m_queue);
  auto in_node_coord = ax::viewIn(command, m_node_coord);
  m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3(cell_lid, cn_cv, in_node_coord); });
}
```

</details>

### Linear Operator Assembly 

The `setValue(DoFLocalId row, DoFLocalId col, Real value)` method of `BSRMatrix` can be used to set coeffcients in the matrix.
You can access the representation of `BSRFormat` i.e. `BSRMatrix` by using the `matrix()` method.

Exemple (applying Dirichlet boundary condition):
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

### Translate BSRFormat to DoFLinearSystem

ArcaneFEM uses `DoFLinearSystem` class to represent and handle linear systems. You can translate the BSRFormat to a given linear system using the `toLinearSystem(DoFLinearSystem &ls)` method:

</details>

<details open>
  <summary><h2>Implementation Details</h2></summary>

  This part is dedicated to the implementation details of `BSRFormat`, in particular the bilinear assembly algorithms.

  ### Default Bilinear Operator Assembly

  #### Construction of the Sparsity

  ##### Building `row_index`

  1. Compute the `neighbors` array: At index `i`, `neighbors` contains the number of neighbors i.e. the number of connected nodes of node `i`.

     a. Iterate over the elements of the mesh in parallel. Store each edge of the element in the `edge` array. Edge `i` (`i: 0 -> nb_edge_per_element`) is stored at index `cur_element_idx * nb_edge_per_element + i`. Edges are represented using a 64-bit integer. The first 32 bits store the `id` of the smaller node in the edge, while the last 32 bits store the `id` of the larger node.

     b. Sort the `edges` array into `sorted_edges`.

     c. Iterate over the `sorted_edges` array in parallel. For each edge, if `sorted_edges[cur_edge_idx + 1] != cur_edge`, increment `neighbors[src]` and `neighbors[dst]` by `1` with an `atomic` operation. This conditional is needed to ensure that we don't count the same edge multiple time (for edges that are shared between multiple elements of the mesh).

  2. `row_index` is the [exclusive scan](https://en.wikipedia.org/wiki/Prefix_sum#Inclusive_and_exclusive_scans) of `neighbors`.

 ##### Building `columns`

 1. Compute `sorted_edges` (see 1.b)
 2. Iterate over the edges in parallel. For each edge, if `sorted_edges[cur_edge_idx + 1] != cur_edge`, "register" the link `src -- dst` in `columns`:

     a. Get the start position of the row `src` in the matrix with `row_index[src]`.

     b. Get the current offset in the `src` row wich `offsets[src]`.

     c. Put `dst` at `columns[start + offset]`.

     d. Increment `offsets[src]` by `1`.

### Atomic-free Bilinear Operator Assembly

#### Construction of the Sparsity

The only step in the precedent algorithm where an atomic operation is used is at 1.c i.e. for computing the `neighbors` array.

To build the sparsity without atomic, we compute the `neighbors` array using [Arcane's node-node connectivity](https://github.com/arcaneframework/framework/pull/1614). This connectivity is computed by Arcane on-demand. It uses node-edge connectivity under the hood and is not accelerated i.e. it will not use GPU to compute the connectivity.

The rest of the algorithm doesn't change apart from the fact that iterations over the edges are done using cell-edge in 3D (and cell-face in 2D) connectivities of Arcane.
