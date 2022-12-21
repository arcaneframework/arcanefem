#ifndef PASSMO_ELASTODYNAMICMODULE_H
#define PASSMO_ELASTODYNAMICMODULE_H

#include "TypesElastodynamic.h"
#include "Elastodynamic_axl.h"
#include "FemUtils.h"
#include "DoFLinearSystem.h"
#include "utilFEM.h"


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::FemUtils;

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
    CellFEMDispatcher m_cell_fem_dispatch;
    Integer m_nb_neqs{0}; // This will be the size of the linear system to solve

private:

    void _initDofs();

    void _applyInputMotion();

    void _computeRHS();

    void _doSolve();

    void _initBoundaryConditions();

    void _initInitialConditions();

//    FixedMatrix<2, 3> _computeBMatrixT3(Cell cell);
    void _applyBoundaryConditions();

    /**
     *  Predict nodal dofs vector (d,v,a) for the Newmark or Generalized-alfa time integration schemes
     */
    void _predictNewmark();

    /**
     *  Update nodal dofs vector (d,v,a) for the Newmark or Generalized-alfa time integration schemes
     */
    void _updateNewmark();
};
#endif // PASSMO_ELASTODYNAMICMODULE_H

