#ifndef PASSMO_ELASTODYNAMICMODULE_H
#define PASSMO_ELASTODYNAMICMODULE_H

#include "TypesElastodynamic.h"
#include "Elastodynamic_axl.h"
#include "FemUtils.h"
#include "utilFEM.h"
#include "DoFLinearSystem.h"
#include "FemDoFsOnNodes.h"


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;
Real REL_PREC{1.0e-15};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Class to define the seismic inputmotion features
 */
/*class Inputmotion{
 public:
  Real3	m_ampli_factors{1.,1.,1.};// Amplification factors to apply on (X, Y, Z) components (default = 1. => no amplification)
  Real  m_max_frequency{10.};// Max frequency for the input signal
  CaseTable* m_acc{ nullptr};
  CaseTable* m_vel{ nullptr};
  CaseTable* m_displ{ nullptr};
  bool m_rigid_base{true};
  bool m_is_vel{ false};
  bool m_is_displ{false};
  bool m_component[3]{true,true,true};//= true when the "ith" component has an input motion curve set
  NodeGroup m_node_group{};

  Inputmotion() = default;
  ~Inputmotion() = default;
};*/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Module Elastodynamic.
 */

class ElastodynamicModule
: public ArcaneElastodynamicObject {
 public:
  explicit ElastodynamicModule(const ModuleBuildInfo &mbi);

 public:

  //! Method called at the beginning of the simulation
  void startInit() override;

  //! Method called at each step
  void compute() override;

  VersionInfo versionInfo() const override;

 private:

  DoFLinearSystem m_linear_system;
  FemDoFsOnNodes m_dofs_on_nodes;

  // Struct to make sure we are using a CaseTable associated
  // to the right file
  struct CaseTableInfo
  {
    String file_name;
    CaseTable* case_table = nullptr;
  };
  // List of CaseTable for traction boundary conditions
  UniqueArray<CaseTableInfo> m_traction_case_table_list;

  // List of CaseTable for dirichlet boundary (surface) conditions
  UniqueArray<CaseTableInfo> m_sacc_case_table_list;
  UniqueArray<CaseTableInfo> m_sdispl_case_table_list;
  UniqueArray<CaseTableInfo> m_svel_case_table_list;
  UniqueArray<CaseTableInfo> m_sforce_case_table_list;

  // List of CaseTable for dirichlet point conditions
  UniqueArray<CaseTableInfo> m_acc_case_table_list;
  UniqueArray<CaseTableInfo> m_displ_case_table_list;
  UniqueArray<CaseTableInfo> m_vel_case_table_list;
  UniqueArray<CaseTableInfo> m_force_case_table_list;

  // List of CaseTable for paraxial boundary conditions
  // (if incident transient wave fields are defined)
  // TO DO ***   UniqueArray<CaseTableInfo> m_paraxial_case_table_list;

  Integer3 integ_order{2, 2, 2};
  Int32 NDIM{2};
  CellFEMDispatcher cell_fem{};
  GaussPointDispatcher gausspt{};
  Real3 gravity{0.,0.,0.};
  Real penalty{1.e64};
  Real gamma{0.5};
  Real beta{0.25};
  Real alfam{0.};
  Real alfaf{0.};
  bool is_alfa_method{false},keep_constop{false};
  Real dt2{0.};
  Int32 linop_nstep{1000}, linop_nstep_counter{0};
  TypesElastodynamic::eElastType elast_type{TypesElastodynamic::NoElastPropType};
  TypesElastodynamic::eAnalysisType analysis_type{TypesElastodynamic::PlaneStrain};

 private:

  void _initDofs();
  void _initCells();
  void _applyInitialNodeConditions();
  void _applyInitialCellConditions();
  // void _applyInputMotion();
  void _assembleLinearLHS();
  void _assembleLinearRHS();
  void _doSolve();
  void _initBoundaryConditions();
  void _applyDirichletBoundaryConditions();
  void _getParaxialContribution(VariableDoFReal& rhs_values);
  void _assembleLHSParaxialContribution();
  void _getTractionContribution(Arcane::VariableDoFReal& rhs_values);
  void _applyNeumannBoundaryConditions();
  Real3x3 _computeJacobian(const ItemWithNodes& cell, const Int32& ig, const RealUniqueArray& vec, Real& jac);
  Real _computeFacLengthOrArea(const Face& face);

  /*  Update nodal dofs vector for the Newmark or Generalized-alfa time integration schemes */
  void _updateNewmark();

  void _computeK(const Cell& cell,const Int32& ig, const RealUniqueArray& vec, const Real3x3& jac, RealUniqueArray2& Ke);
  void _computeElemMass(const Cell& cell,const Int32& ig, const RealUniqueArray& vec, const Real& jacobian, RealUniqueArray2& Me);
  void _computeKParax(const Face& face, const Int32& ig, const RealUniqueArray& vec, const Real& jacobian,
                      RealUniqueArray2& Ke, const Real3& RhoC);
  void _computeTracParax(const Face& face, const Int32& ig, const RealUniqueArray& vec, const Real& jacobian,
                         RealUniqueArray& Fe, const Real3& RhoC);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif // PASSMO_ELASTODYNAMICMODULE_H
