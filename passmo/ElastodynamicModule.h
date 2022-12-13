#ifndef PASSMO_ELASTODYNAMICMODULE_H
#define PASSMO_ELASTODYNAMICMODULE_H

#include "TypesElastodynamic.h"
#include "Elastodynamic_axl.h"
#include "FemUtils.h"
#include "FemLinearSystemDoF.h"


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
using namespace Arcane;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Module Elastodynamic.
 */

class ElastodynamicModule
: public ArcaneElastodynamicObject
    {
    public:
        explicit ElastodynamicModule(const ModuleBuildInfo& mbi);

    public:

    //! Method called at the beginning of the simulation
    void startInit() override;

    //! Method called at each step
    void predict() override;

    //! Method called at each step
    void compute() override;

    //! Method called at each step
    void update() override;

    VersionInfo versionInfo() const override;

private:
    FemLinearSystemDoF m_linear_system;

private:

    void _applyInputMotion();
    void _computeRHS();
    void _doSolve();
    void _initBoundaryConditions();
    void _initInitialConditions();
//    FixedMatrix<2, 3> _computeBMatrixT3(Cell cell);
    void _applyBoundaryConditions();
    };

#endif // PASSMO_ELASTODYNAMICMODULE_H

