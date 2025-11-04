// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* FEM code to test vectorial FE for Elastodynamics problem.                 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"
#include "Traction.h"
#include "Dirichlet.h"
#include "SourceTerm.h"

/*---------------------------------------------------------------------------*/
/**
 * @brief Initializes the FemModule at the start of the simulation.
 *
 * This method initializes degrees of freedom (DoFs) on nodes.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "[ArcaneFem-Info] Started module  startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_nodes.initialize(mesh(), mesh()->dimension());

  _getParameters();

  bool use_csr_in_linearsystem = options()->linearSystem.serviceName() == "HypreLinearSystem";
  if (m_matrix_format == "BSR")
    m_bsr_format.initialize(defaultMesh(), mesh()->dimension(), use_csr_in_linearsystem, 0);
  else if (m_matrix_format == "AF-BSR")
    m_bsr_format.initialize(defaultMesh(), mesh()->dimension(), use_csr_in_linearsystem, 1);

  t = dt;
  tmax = tmax - dt;
  m_global_deltat.assign(dt);

  _readCaseTables();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*
 * @brief Computes the FEM simulation for all time step.
 *
 * This method performs the main computation of the FEM simulation.
 * It assembles the linear system, solves it, and updates the variables.
 * It also handles stopping the computation loop if the maximum time is reached.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "[ArcaneFem-Info] Started module  compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop the computation loop if the maximum time is reached
  if (t >= tmax)
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  info() << "Time iteration at t : " << t << " (s) ";

  // Set if we want to keep the matrix structure between calls.
  // The matrix has to have the same structure (same structure for non-zero)
  bool keep_struct = true;
  if (m_linear_system.isInitialized() && keep_struct) {
    m_linear_system.clearValues();
  }
  else {
    m_linear_system.reset();
    m_linear_system.setLinearSystemFactory(options()->linearSystem());
    m_linear_system.initialize(subDomain(),  acceleratorMng()->defaultRunner(), m_dofs_on_nodes.dofFamily(), "Solver");
  }

  if (m_petsc_flags != NULL) {
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();
  _updateTime();

  if ((t > tmax + dt - 1e-8) && m_cross_validation)
    _validateResults();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateTime()
{
  t += dt;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_doStationarySolve()
{
  if(m_assemble_linear_system){
    _assembleBilinearOperator();
    _assembleLinearOperator();
  }
  if(m_solve_linear_system){
    _solve();
    _updateVariables();
  }
}

/*---------------------------------------------------------------------------*/
/*
 * @brief Reads the parameters from the options.
 *
 * This method retrieves various parameters such as time step, material properties,
 * and time discretization methods from the options provided in the simulation.
 * It also initializes the material properties based on the retrieved parameters.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_getParameters()
{
  info() << "[ArcaneFem-Info] Started module  _getParameters()";
  Real elapsedTime = platform::getRealTime();

  //--------- time parameters -----------//
  tmax = options()->tmax(); // max time 𝑡ₘₐₓ
  dt = options()->dt(); // time step δ𝑡

  //--- damping term parameter ---//
  etam = options()->etam(); // damping parameter ηₘ
  etak = options()->etak(); // damping parameter ηₖ

  //--- time discretization parameter ---//
  alpm = options()->alpm(); // time discretization parameter αᵐ
  alpf = options()->alpf(); // time discretization parameter ηᶠ

  //--------- material parameter ---------//
  E = options()->E(); // Youngs modulus 𝐸
  nu = options()->nu(); // Poisson ratio ν
  rho = options()->rho(); // Density ρ

  //--------- mesh parameter ------------//
  m_hex_quad_mesh = options()->hexQuadMesh();

  mu = E / (2 * (1 + nu)); // lame parameter μ
  lambda = E * nu / ((1 + nu) * (1 - 2 * nu)); // lame parameter λ

  if( options()->mu.isPresent())
    mu = options()->mu;

  if( options()->lambda.isPresent())
    lambda = options()->lambda;

  //--------- body force ---------//
  if (options()->f.isPresent()) {
    const UniqueArray<String> f_string = options()->f();
    info() << "[ArcaneFem-Info] Applying Bodyforce " << f_string;
    for (Int32 i = 0; i < f_string.size(); ++i) {
      if (f_string[i] != "NULL") {
        f[i] = std::stod(f_string[i].localstr());
      }
    }
  }

  //----- time discretization Newmark-Beta or Generalized-alpha  -----//
  if (options()->timeDiscretization == "Newmark-beta") {

    info() << "[ArcaneFem-Info] Apply time discretization via Newmark-beta ";

    gamma = 0.5;
    beta  = (1./4.)*(gamma+0.5)*(gamma+0.5)  ;

    c0 = rho / (beta * dt * dt) + etam * rho * gamma / beta / dt;
    c1 = lambda + lambda * etak * gamma / beta / dt;
    c2 = mu + mu * etak * gamma / beta / dt;
    c3 = rho / beta / dt - etam * rho * (1 - gamma / beta);
    c4 = rho * ((1. - 2. * beta) / 2. / beta - etam * dt * (1. - gamma / 2 / beta));
    c5 = -lambda * etak * gamma / beta / dt;
    c6 = -mu * etak * gamma / beta / dt;
    c7 = etak * lambda * (gamma / beta - 1);
    c8 = etak * lambda * dt * ((1. - 2 * beta) / 2. / beta - (1. - gamma));
    c9 = etak * mu * (gamma / beta - 1);
    c10 = etak * mu * dt * ((1. - 2 * beta) / 2. / beta - (1. - gamma));
  }
  else if (options()->timeDiscretization == "Generalized-alpha") {

    info() << "[ArcaneFem-Info] Apply time discretization via Generalized-alpha ";

    gamma = 0.5 + alpf - alpm                ;
    beta  = (1./4.)*(gamma+0.5)*(gamma+0.5)  ;

    c0 = rho * (1. - alpm) / (beta * dt * dt) + etam * rho * gamma * (1 - alpf) / beta / dt;
    c1 = lambda * (1. - alpf) + lambda * etak * gamma * (1. - alpf) / beta / dt;
    c2 = mu * (1. - alpf) + mu * etak * gamma * (1. - alpf) / beta / dt;
    c3 = rho * (1. - alpm) / beta / dt - etam * rho * (1 - gamma * (1 - alpf) / beta);
    c4 = rho * ((1. - alpm) * (1. - 2. * beta) / 2. / beta - alpm - etam * dt * (1. - alpf) * (1. - gamma / 2 / beta));
    c5 = lambda * alpf - lambda * etak * gamma * (1. - alpf) / beta / dt;
    c6 = mu * alpf - mu * etak * gamma * (1. - alpf) / beta / dt;
    c7 = etak * lambda * (gamma * (1. - alpf) / beta - 1);
    c8 = etak * lambda * dt * (1. - alpf) * ((1. - 2 * beta) / 2. / beta - (1. - gamma));
    c9 = etak * mu * (gamma * (1. - alpf) / beta - 1);
    c10 = etak * mu * dt * (1. - alpf) * ((1. - 2 * beta) / 2. / beta - (1. - gamma));
  }
  else {
    ARCANE_FATAL("Only Newmark-beta | Generalized-alpha are supported for time-discretization ");
  }

  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->crossValidation();
  m_petsc_flags = options()->petscFlags();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*
  * @brief Reads case tables for traction boundary conditions.
  *
  * This method reads the case tables specified in the options and stores
  * them in a list for later use.
  */
/*---------------------------------------------------------------------------*/

void FemModule::
_readCaseTables()
{
  IParallelMng* pm = subDomain()->parallelMng();
  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  // loop over all traction boundries
  for (BC::ITractionBoundaryCondition* bs : bc->tractionBoundaryConditions()) {
    CaseTable* case_table = nullptr;
    auto traction_table_file_name = bs->getTractionInputFile();
    bool getTractionFromTable = !traction_table_file_name.empty();
    if (getTractionFromTable)
      case_table = readFileAsCaseTable(pm, traction_table_file_name, 3);
    m_traction_case_table_list.add(CaseTableInfo{ traction_table_file_name, case_table });
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateVariables()
{
  info() << "[ArcaneFem-Info] Started module  _updateVariables()";
  Real elapsedTime = platform::getRealTime();

  VariableDoFReal& dof_u(m_linear_system.solutionVariable());
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // First loop: Update displacement values
  const Int8 dim = mesh()->dimension();
  ENUMERATE_ (Node, inode, ownNodes()) {
    Node node = *inode;
    Real3 u_disp = Real3::zero();
    for (Int8 i = 0; i < dim; ++i) {
      u_disp[i] = dof_u[node_dof.dofId(node, i)];
    }
    m_dU[node] = u_disp;
  }

  // Synchronize the variables
  m_dU.synchronize();
  m_U.synchronize();
  m_V.synchronize();
  m_A.synchronize();

  // Second loop: Update acceleration, velocity and displacement
  ENUMERATE_ (Node, inode, allNodes()) {
    Node node = *inode;
    Real3 aloc = Real3::zero(); // Local acceleration vector

    // Constants for Newmark-beta method
    const Real beta_dt2 = beta * (dt * dt);
    const Real beta_factor = (1. - 2. * beta) / (2. * beta);
    const Real gamma_factor = (1. - gamma);

    // Update each component
    for (Integer i = 0; i < dim; ++i) {
      // Calculate new acceleration
      aloc[i] = (m_dU[node][i] - m_U[node][i] - dt * m_V[node][i]) / beta_dt2 - beta_factor * m_A[node][i];

      // Update velocity, acceleration, and displacement
      m_V[node][i] += dt * (gamma_factor * m_A[node][i] + gamma * aloc[i]);
      m_A[node][i] = aloc[i];
      m_U[node][i] = m_dU[node][i];
    }
  }

    elapsedTime = platform::getRealTime() - elapsedTime;
    ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "update-variables", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*
 * @brief Assembles the linear operator for the FEM simulation.
 *
 * This method computes the linear operator for the FEM simulation
 * by assembling the linear parts of the system. It uses the
 * material properties and the mesh to compute the contributions
 * from each element. The results are stored in the right-hand side
 * of the linear system.
 * 
 * @note This method assumes that the mesh and material properties
 *  are already set up.
 */
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  //if(t <= dt) // this works for PETSc but breaks Hypre
  if(m_matrix_format == "BSR" || m_matrix_format == "AF-BSR")
    m_bsr_format.toLinearSystem(m_linear_system);

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  _applySourceTerm(rhs_values, node_dof);
  _applyTraction(rhs_values, node_dof);
  _applyDirichlet(rhs_values, node_dof);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "rhs-vector-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTria3Gpu()
{
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto command = makeCommand(acceleratorMng()->defaultQueue());
  auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
  auto c0_copy = c0;
  auto c1_copy = c1;
  auto c2_copy = c2;

  m_bsr_format.computeSparsity();
  if (m_matrix_format == "BSR")
    m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTria3Gpu(cell_lid, cn_cv, in_node_coord, c0_copy, c1_copy, c2_copy); });
  else
    m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTria3Gpu(cell_lid, cn_cv, in_node_coord, c0_copy, c1_copy, c2_copy, node_lid); });
  }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperatorTetra4Gpu()
{
  UnstructuredMeshConnectivityView m_connectivity_view(mesh());
  auto cn_cv = m_connectivity_view.cellNode();
  auto command = makeCommand(acceleratorMng()->defaultQueue());
  auto in_node_coord = Accelerator::viewIn(command, m_node_coord);
  auto c0_copy = c0;
  auto c1_copy = c1;
  auto c2_copy = c2;

  m_bsr_format.computeSparsity();
  if (m_matrix_format == "BSR")
    m_bsr_format.assembleBilinearAtomic([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid) { return computeElementMatrixTetra4Gpu(cell_lid, cn_cv, in_node_coord, c0_copy, c1_copy, c2_copy); });
  else
    m_bsr_format.assembleBilinearAtomicFree([=] ARCCORE_HOST_DEVICE(CellLocalId cell_lid, Int32 node_lid) { return computeElementVectorTetra4Gpu(cell_lid, cn_cv, in_node_coord, c0_copy, c1_copy, c2_copy, node_lid); });
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Assembles the FEM bilinear operator on CPU.
 *
 * This method assembles the FEM stiffness matrix by iterating over each cell,
 * computing the element stiffness matrix using the provided function, and
 * populating the global stiffness matrix accordingly.
 *
 * @tparam N Total DOF size (nodes_per_element × dimensions).
 * @param compute_element_matrix function computing cell's element stiffness matrix.
 */
/*---------------------------------------------------------------------------*/

template <int N>
void FemModule::
_assembleBilinearOperatorCpu(const std::function<RealMatrix<N, N>(const Cell&)>& compute_element_matrix)
{
  const Int32 dim = mesh()->dimension();
  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  ENUMERATE_ (Cell, icell, allCells()) {
    Cell cell = *icell;
    auto K_e = compute_element_matrix(cell);

    Int32 n1_index = 0;
    for (Node node1 : cell.nodes()) {
      if (node1.isOwn()) {
        Int32 n2_index = 0;
        for (Node node2 : cell.nodes()) {
          for (Int32 i = 0; i < dim; ++i) {
            DoFLocalId dof1 = node_dof.dofId(node1, i);
            for (Int32 j = 0; j < dim; ++j) {
              DoFLocalId dof2 = node_dof.dofId(node2, j);
              Real value = K_e(dim * n1_index + i, dim * n2_index + j);
              m_linear_system.matrixAddValue(dof1, dof2, value);
            }
          }
          ++n2_index;
        }
      }
      ++n1_index;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleBilinearOperator()
{
  info() << "[ArcaneFem-Info] Started module  _assembleBilinearOperator()";
  Real elapsedTime = platform::getRealTime();

  if (t <= dt) {
    if (m_matrix_format == "DOK") {
      if (mesh()->dimension() == 2) {
        if (m_hex_quad_mesh) {
          _assembleBilinearOperatorCpu<8>([this](const Cell& cell) { return _computeElementMatrixQuad4(cell); });
        }
        else {
          _assembleBilinearOperatorCpu<6>([this](const Cell& cell) { return _computeElementMatrixTria3(cell); });
        }
      }
      if (mesh()->dimension() == 3) {
        if (m_hex_quad_mesh) {
          _assembleBilinearOperatorCpu<24>([this](const Cell& cell) { return _computeElementMatrixHexa8(cell); });
        }
        else{
          _assembleBilinearOperatorCpu<12>([this](const Cell& cell) { return _computeElementMatrixTetra4(cell); });
        }
      }
    }
    else if (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR") {
      if (mesh()->dimension() == 2) {
        _assembleBilinearOperatorTria3Gpu();
      }
      if (mesh()->dimension() == 3) {
        _assembleBilinearOperatorTetra4Gpu();
      }
    }
    else {
      ARCANE_FATAL("Unsupported matrix type, only DOK| BSR | AF-BSR is supported.");
    }
  }

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "lhs-matrix-assembly", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_solve()
{
  info() << "[ArcaneFem-Info] Started module  _solve()";
  Real elapsedTime = platform::getRealTime();

  m_linear_system.applyLinearSystemTransformationAndSolve();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "solve-linear-system", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_validateResults()
{
  info() << "[ArcaneFem-Info] Started module  _validateResults()";
  Real elapsedTime = platform::getRealTime();

  if (allNodes().size() < 200) {
    int p = std::cout.precision();
    std::cout.precision(17);
    ENUMERATE_ (Node, inode, allNodes()) {
      Node node = *inode;
      std::cout << "( N_id, u1, u2, u3 ) = ( "
                << node.uniqueId() << ", " << m_dU[node].x << ", " << m_dU[node].y << ", " << m_dU[node].z
                << ")\n";
    }
    std::cout.precision(p);
  }

  String filename = options()->resultFile();
  const double epsilon = 1.0e-4;
  const double min_value_to_test = 1.0e-14;

  info() << "[ArcaneFem-Info] Validating results filename=" << filename << " epsilon =" << epsilon;

  if (!filename.empty())
    Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_dU, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "cross-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
