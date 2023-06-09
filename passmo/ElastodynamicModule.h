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

   // List of CaseTable for paraxial boundary conditions
   // (if incident transient wave fields are defined)
// TO DO ***   UniqueArray<CaseTableInfo> m_paraxial_case_table_list;

   Integer3 integ_order{2, 2, 2};
   Int32 NDIM{2};
   CellFEMDispatcher cell_fem{};
   GaussPointDispatcher gausspt{};
   Real3 gravity{0.,0.,-9.81};
   Real penalty{1.e30};
   Real gamma{0.5};
   Real beta{0.25};
   Real alfam{0.};
   Real alfaf{0.};
   bool is_alfa_method{false},keep_constop{false};
   Real dt2{0.};
   Int32 linop_nstep{100}, linop_nstep_counter{0};

private:

 void _initDofs();
 void _applyInitialNodeConditions();
 void _applyInitialCellConditions();
// void _applyInputMotion();
 void _assembleLinearGlobal2D();
 void _assembleLinearGlobal3D();
 void _doSolve();
 void _initBoundaryConditions();
 //void _initInputMotion();
 void _applyBoundaryConditions();
 Real3x3 _computeInverseJacobian3D(const Cell& cell,const Real3& ref_coord, Real& jacobian);
 Real2x2  _computeInverseJacobian2D(const Cell& cell,const Real3& ref_coord, Real& jacobian);
 Real _computeFacLengthOrArea(const Face& face);


 //    FixedMatrix<2, 3> _computeBMatrixT3(Cell cell);
 /**
     *  Predict nodal dofs vector for the Newmark or Generalized-alfa time integration schemes
     */
 void _predictNewmark();
 /**
     *  Update nodal dofs vector for the Newmark or Generalized-alfa time integration schemes
     */
 void _updateNewmark();

 void _computeKMF3D(const Cell& cell,RealUniqueArray2& Ke, RealUniqueArray2& Me, RealUniqueArray& Fe);
 void _computeKMF2D(const Cell& cell,RealUniqueArray2& Ke, RealUniqueArray2& Me, RealUniqueArray& Fe);
};

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/

#endif // PASSMO_ELASTODYNAMICMODULE_H

