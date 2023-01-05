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

private:

 void _initDofs();
 void _applyInitialNodeConditions();
 void _applyInitialCellConditions();
 void _applyInputMotion();
 void _assembleLinearGlobal();
 void _doSolve();
 void _initBoundaryConditions();
 void _applyBoundaryConditions();
 static Real _computeJacobian(const Real3x3& jacmat, const Integer& ndim);
 Real3x3 _computeJacobianMatrix(const Cell& cell,const Real3& ref_coord);
 static Real3x3 _computeInverseJacobianMatrix(const Real3x3& jacmat, const Real& jacobian, const Integer& ndim);

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


#endif // PASSMO_ELASTODYNAMICMODULE_H

