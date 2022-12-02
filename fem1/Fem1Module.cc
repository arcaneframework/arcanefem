// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "Fem1_axl.h"
#include <arcane/ITimeLoopMng.h>

using namespace Arcane;

/*!
 * \brief Module Fem1.
 */
class Fem1Module
: public ArcaneFem1Object
{
 public:
  explicit Fem1Module(const ModuleBuildInfo& mbi) 
  : ArcaneFem1Object(mbi) { }

 public:
  /*!
   * \brief Méthode appelée à chaque itération.
   */
  void compute() override;
  /*!
   * \brief Méthode appelée lors de l'initialisation.
   */
  void startInit() override;

  /** Retourne le numéro de version du module */
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 private:

  void _doStationarySolve();
  void _updateBoundayConditions();
  void _computeConductivity();
  void _computeGeneralizedFluxes();
  void _solve();
  void _initBoundaryconditions();
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
compute()
{
  info() << "Module Fem1 COMPUTE";

  // Stop code after computatio
  subDomain()->timeLoopMng()->stopComputeLoop(true);

  info() << "NB_CELL=" << allCells().size() << " NB_FACE=" << allFaces().size();
  _doStationarySolve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
startInit()
{
  info() << "Module Fem1 INIT";

  // # init mesh
  // # init behavior
  // # init behavior on mesh entities
  // # init BCs
  _initBoundaryconditions();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_doStationarySolve()
{
  // # update BCs
  _updateBoundayConditions();

  // K=self.compute_conductivity()
  _computeConductivity();
  //      # compute RHS
  _computeGeneralizedFluxes();

  // # adapt K and RHS to take into account Dirichlet BCs
  //         for node in self.mesh.nodes:
  //             if node.is_T_fixed:
  //                 K[node.rank,node.rank]=K[node.rank,node.rank]+10**6
  //                 RHS[node.rank]=RHS[node.rank]+(10**6)*node.T

  info() << "TODO: adapt K and RHS";

  // # T=linalg.solve(K,RHS)
  _solve();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_initBoundaryconditions()
{
  info() << "Init boundary conditions...";

  info() << "TODO: Init boundary conditions";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_updateBoundayConditions()
{
  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_computeConductivity()
{
  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_computeGeneralizedFluxes()
{
  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fem1Module::
_solve()
{
  info() << "TODO " << A_FUNCINFO;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_FEM1(Fem1Module);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
