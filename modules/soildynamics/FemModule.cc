// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2025 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* FemModule.cc                                                (C) 2022-2025 */
/*                                                                           */
/* FEM code to test vectorial FE for Soildynamics problem.                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "FemModule.h"
#include "ElementMatrix.h"
#include "ElementMatrixHexQuad.h"
#include "SourceTerm.h"
#include "Dirichlet.h"
#include "Traction.h"
#include "Paraxial.h"
#include "DoubleCouple.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
compute()
{
  info() << "[ArcaneFem-Info] Started module  compute()";
  Real elapsedTime = platform::getRealTime();

  // Stop code after computations
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

  if (m_petsc_flags != NULL){
    CommandLineArguments args = ArcaneFemFunctions::GeneralFunctions::getPetscFlagsFromCommandline(m_petsc_flags);
    m_linear_system.setSolverCommandLineArguments(args);
  }

  _doStationarySolve();
  if (t >= tmax  && m_cross_validation)
    _validateResults();
  _updateTime();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "compute", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
startInit()
{
  info() << "[ArcaneFem-Info] Started module  startInit()";
  Real elapsedTime = platform::getRealTime();

  m_dofs_on_nodes.initialize(mesh(), mesh()->dimension());

  _getParameters();

  bool use_csr_in_linear_system =
    options()->linearSystem.serviceName() == "HypreLinearSystem" ||
    options()->linearSystem.serviceName() == "PETScLinearSystem";

  if (m_matrix_format == "BSR")
    m_bsr_format.initialize(defaultMesh(), mesh()->dimension(), use_csr_in_linear_system, 0);
  else if (m_matrix_format == "AF-BSR")
    m_bsr_format.initialize(defaultMesh(), mesh()->dimension(), use_csr_in_linear_system, 1);

  t = dt;
  tmax = tmax;
  m_global_deltat.assign(dt);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "initialize", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_updateTime()
{
  t += dt;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Performs a stationary solve for the FEM system.
 *
 * This method follows via the following steps:
 *   1. _assembleBilinearOperator()  Assembles the FEM  matrix 𝐀
 *   2. _assembleLinearOperator()    Assembles the FEM RHS vector 𝐛
 *   3. _solve()                     Solves for solution vector 𝐮 = 𝐀⁻¹𝐛
 *   4. _updateVariables()           Updates FEM variables 𝐮 = 𝐱
 */
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
/*---------------------------------------------------------------------------*/

void FemModule::
_getParameters()
{
  info() << "[ArcaneFem-Info] Started module  _getParameters()";
  Real elapsedTime = platform::getRealTime();

  //--------- mesh parameter ------------//
  m_hex_quad_mesh = options()->hexQuadMesh();

  //--------- time parameters -----------//
  tmax = options()->tmax(); // max time
  dt = options()->dt(); // time step

  //--------- material parameter ---------//
  E = options()->E(); // Youngs modulus
  nu = options()->nu(); // Poission ratio
  rho = options()->rho(); // Density
  cp = options()->cp(); // Wave velocity primary
  cs = options()->cs(); // Wave velocity secondary

  if (options()->E.isPresent() && options()->nu.isPresent()) {
    mu = E / (2 * (1 + nu));
    lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    cs = math::sqrt(mu / rho);
    cp = math::sqrt((lambda + (2. * mu)) / rho);
  }

  if (options()->mu.isPresent() && options()->lambda.isPresent()) {
    mu = options()->mu;
    lambda = options()->lambda;
    cs = math::sqrt(mu / rho);
    cp = math::sqrt((lambda + (2. * mu)) / rho);
  }

  if ((options()->cp.isPresent()) && (options()->cs.isPresent())) {
    mu = cs * cs * rho;
    lambda = cp * cp * rho - 2 * mu;
  }

  if (options()->f.isPresent()) {
    const UniqueArray<String> f_string = options()->f();
    info() << "[ArcaneFem-Info] Applying Bodyforce " << f_string;
    for (Int32 i = 0; i < f_string.size(); ++i) {
      if (f_string[i] != "NULL") {
        f[i] = std::stod(f_string[i].localstr());
      }
    }
  }

  gamma = 0.5;
  beta = (1. / 4.) * (gamma + 0.5) * (gamma + 0.5);

  c0 = rho / (beta * dt * dt);
  c1 = lambda;
  c2 = mu;
  c3 = rho / (beta * dt);
  c4 = rho * (1. / 2. / beta - 1.);
  c5 = 0.;
  c6 = 0.;
  c7 = rho * gamma / beta / dt;
  c8 = rho * (1. - gamma / beta);
  c9 = rho * dt * (1. - gamma / (2. * beta));

  m_matrix_format = options()->matrixFormat();
  m_assemble_linear_system = options()->assembleLinearSystem();
  m_solve_linear_system = options()->solveLinearSystem();
  m_cross_validation = options()->hasSolutionComparisonFile();
  m_petsc_flags = options()->petscFlags();

  _readCaseTables();

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "get-material-params", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void FemModule::
_readCaseTables()
{
  IParallelMng* pm = subDomain()->parallelMng();
  BC::IArcaneFemBC* bc = options()->boundaryConditions();

  // loop over all traction boundries
  if (bc) {
    for (BC::ITractionBoundaryCondition* bs : bc->tractionBoundaryConditions()) {
      CaseTable* case_table = nullptr;
      auto traction_table_file_name = bs->getTractionInputFile();
      bool getTractionFromTable = !traction_table_file_name.empty();
      if (getTractionFromTable)
        case_table = readFileAsCaseTable(pm, traction_table_file_name, 3);
      m_traction_case_table_list.add(CaseTableInfo{ traction_table_file_name, case_table });
    }
  }

  // loop over all double couple boundries
  for (const auto& bs : options()->doubleCouple()) {
    CaseTable* case_table = nullptr;
    if (bs->doubleCoupleInputFile.isPresent())
      case_table = readFileAsCaseTable(pm, bs->doubleCoupleInputFile(), 1);
    m_double_couple_case_table_list.add(Arcane::FemUtils::CaseTableInfo{ bs->doubleCoupleInputFile(), case_table });
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
/*---------------------------------------------------------------------------*/

void FemModule::
_assembleLinearOperator()
{
  info() << "[ArcaneFem-Info] Started module _assembleLinearOperator()";
  Real elapsedTime = platform::getRealTime();

  // Determine if using BSR format
  bool use_bsr = (m_matrix_format == "BSR" || m_matrix_format == "AF-BSR");
  BSRMatrix* bsr_matrix = use_bsr ? &(m_bsr_format.matrix()) : nullptr;

  // Temporary variable to keep values for the RHS part of the linear system
  VariableDoFReal& rhs_values(m_linear_system.rhsVariable());
  rhs_values.fill(0.0);

  auto node_dof(m_dofs_on_nodes.nodeDoFConnectivityView());

  // Apply all boundary conditions and source terms
  _applySourceTerm(rhs_values, node_dof);
  _applyTraction(rhs_values, node_dof);
  _applyParaxial(rhs_values, node_dof, bsr_matrix);
  _applyDoubleCouple(rhs_values, node_dof);
  _applyDirichlet(rhs_values, node_dof);

  // Convert BSR format to linear system if needed
  if (use_bsr)
    m_bsr_format.toLinearSystem(m_linear_system);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "rhs-vector-assembly", elapsedTime);
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

  String filename = options()->solutionComparisonFile();
  const double epsilon = 1.0e-4;
  const double min_value_to_test = 1.0e-10;

  Arcane::FemUtils::checkNodeResultFile(traceMng(), filename, m_dU, epsilon, min_value_to_test);

  elapsedTime = platform::getRealTime() - elapsedTime;
  ArcaneFemFunctions::GeneralFunctions::printArcaneFemTime(traceMng(), "cross-validation", elapsedTime);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM(FemModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
