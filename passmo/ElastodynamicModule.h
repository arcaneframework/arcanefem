#ifndef PASSMO_ELASTODYNAMICMODULE_H
#define PASSMO_ELASTODYNAMICMODULE_H

#include "TypesElastodynamic.h"
#include "Elastodynamic_axl.h"
#include "FemUtils.h"
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
class Inputmotion{
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
};

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
   Inputmotion m_input{};

private:

 void _initDofs();
 void _applyInitialNodeConditions();
 void _applyInitialCellConditions();
 void _applyInputMotion();
 void _assembleLinearGlobal();
 void _doSolve();
 void _initBoundaryConditions();
 void _initInputMotion();
 void _applyBoundaryConditions();
 static Real _computeJacobian(const Real3x3& jacmat, const Integer& ndim);
 Real3x3 _computeJacobianMatrix(const Cell& cell,const Real3& ref_coord);
 static Real3x3 _computeInverseJacobianMatrix(const Real3x3& jacmat, const Real& jacobian, const Integer& ndim);
 static Real _interpolCurve(const Real& y_before, const Real& y_after,
                const Real& x, const Real& x_before, const Real& x_after);

 //    FixedMatrix<2, 3> _computeBMatrixT3(Cell cell);
 /**
     *  Predict nodal dofs vector for the Newmark or Generalized-alfa time integration schemes
     */
 void _predictNewmark();
 /**
     *  Update nodal dofs vector for the Newmark or Generalized-alfa time integration schemes
     */
 void _updateNewmark();

 void _computeKMF(const Cell& cell,RealUniqueArray2& Ke, RealUniqueArray2& Me, RealUniqueArray& Fe);
};

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/

#endif // PASSMO_ELASTODYNAMICMODULE_H

