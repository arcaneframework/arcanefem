#include "arcane/MathUtils.h"
#include <arcane/utils/NumArray.h>
#include <arcane/IParallelMng.h>
#include <arcane/IMesh.h>
#include <arcane/IItemFamily.h>
#include <arcane/geometry/IGeometry.h>
#include <arcane/VariableTypes.h>
#include "IArcaneFemBC.h"
#include "ArcaneFemFunctions.h"

using namespace Arcane;
/*---------------------------------------------------------------------------*/
/**
 * @brief Initialize the CellFEMDispatcher class (coming from PASSMO)
 */
/*---------------------------------------------------------------------------*/
ArcaneFemFunctions::CellFEMDispatcher::CellFEMDispatcher(){
 // Setting to null default value
 for(int i = 0; i < NB_BASIC_ITEM_TYPE; ++i )
 {
   m_shapefunc[i] = nullptr;
   m_shapefuncderiv[i] = nullptr;
 }

 // Gives functions to compute shape function value in finite-element reference coordinate system
 // Linear elements
 m_shapefunc[IT_Line2] = ArcaneFemFunctions::FemShapeMethods::line2ShapeFuncVal;
 m_shapefunc[IT_Triangle3] = ArcaneFemFunctions::FemShapeMethods::tri3ShapeFuncVal;
 m_shapefunc[IT_Quad4] = ArcaneFemFunctions::FemShapeMethods::quad4ShapeFuncVal;
 m_shapefunc[IT_Tetraedron4] = ArcaneFemFunctions::FemShapeMethods::tetra4ShapeFuncVal;
 m_shapefunc[IT_Hexaedron8] = ArcaneFemFunctions::FemShapeMethods::hexa8ShapeFuncVal;
 m_shapefunc[IT_Pentaedron6] = ArcaneFemFunctions::FemShapeMethods::penta6ShapeFuncVal;
 m_shapefunc[IT_Pyramid5] = ArcaneFemFunctions::FemShapeMethods::pyramid5ShapeFuncVal;

 // Quadratic elements
 m_shapefunc[IT_Line3] = ArcaneFemFunctions::FemShapeMethods::line3ShapeFuncVal;
 m_shapefunc[IT_Triangle6] = ArcaneFemFunctions::FemShapeMethods::tri6ShapeFuncVal;
 m_shapefunc[IT_Quad8] = ArcaneFemFunctions::FemShapeMethods::quad8ShapeFuncVal;
 m_shapefunc[IT_Tetraedron10] = ArcaneFemFunctions::FemShapeMethods::tetra10ShapeFuncVal;
 m_shapefunc[IT_Hexaedron20] = ArcaneFemFunctions::FemShapeMethods::hexa20ShapeFuncVal;

 // Gives functions to compute shape function derivate vector at all nodes of a finite-element
 // along a local direction (in reference coordinate system)
 // Linear elements
 m_shapefuncderiv[IT_Line2] = ArcaneFemFunctions::FemShapeMethods::line2ShapeFuncDeriv;
 m_shapefuncderiv[IT_Triangle3] = ArcaneFemFunctions::FemShapeMethods::tri3ShapeFuncDeriv;
 m_shapefuncderiv[IT_Quad4] = ArcaneFemFunctions::FemShapeMethods::quad4ShapeFuncDeriv;
 m_shapefuncderiv[IT_Tetraedron4] = ArcaneFemFunctions::FemShapeMethods::tetra4ShapeFuncDeriv;
 m_shapefuncderiv[IT_Hexaedron8] = ArcaneFemFunctions::FemShapeMethods::hexa8ShapeFuncDeriv;
 m_shapefuncderiv[IT_Pentaedron6] = ArcaneFemFunctions::FemShapeMethods::penta6ShapeFuncDeriv;
 m_shapefuncderiv[IT_Pyramid5] = ArcaneFemFunctions::FemShapeMethods::pyramid5ShapeFuncDeriv;

 // Quadratic elements
 m_shapefuncderiv[IT_Line3] = ArcaneFemFunctions::FemShapeMethods::line3ShapeFuncDeriv;
 m_shapefuncderiv[IT_Triangle6] = ArcaneFemFunctions::FemShapeMethods::tri6ShapeFuncDeriv;
 m_shapefuncderiv[IT_Quad8] = ArcaneFemFunctions::FemShapeMethods::quad8ShapeFuncDeriv;
 m_shapefuncderiv[IT_Tetraedron10] = ArcaneFemFunctions::FemShapeMethods::tetra10ShapeFuncDeriv;
 m_shapefuncderiv[IT_Hexaedron20] = ArcaneFemFunctions::FemShapeMethods::hexa20ShapeFuncDeriv;

}

/*---------------------------------------------------------------------------*/
/**
 * @brief Provides at once, all Gauss data of a given input finite element
 * (vector containing weights, ref. coordinates, nodal shape values & derivatives)
 * This method is generic (coming from PASSMO)
 */
/*---------------------------------------------------------------------------*/
RealUniqueArray ArcaneFemFunctions::CellFEMDispatcher::
getGaussData(ItemWithNodes item, Integer nint, Integer ngauss){

 auto nnod = item.nbNode();
 auto cell_type = item.type();
 ngauss = ArcaneFemFunctions::FemGaussQuadrature::getNbGaussPointsfromOrder(cell_type,nint);

 // Vector of double containing:
 // ngauss points * [weight, gauss ref coord [Real3], nnod * (shapefunc values, 3*shapefunc deriv
 // in ref. coord system)]
 Integer nsize = ngauss * 4 * (1 + nnod);
 RealUniqueArray vec(nsize);

 Integer index{ 0 };
 for (Integer ig = 0; ig < ngauss; ++ig) {
   auto wt = ArcaneFemFunctions::FemGaussQuadrature::getGaussWeight(item, nint, ig);
   auto pos = ArcaneFemFunctions::FemGaussQuadrature::getGaussRefPosition(item, nint, ig);
   vec[index++] = wt;
   vec[index++] = pos.x;
   vec[index++] = pos.y;
   vec[index++] = pos.z;

   for (Int32 inod = 0; inod < nnod; ++inod) {
     auto Phi_i = getShapeFuncVal(cell_type, inod, pos);
     vec[index++] = Phi_i;
     auto dPhi = getShapeFuncDeriv(cell_type, inod, pos);
     vec[index++] = dPhi.x;
     vec[index++] = dPhi.y;
     vec[index++] = dPhi.z;
   }
 }
 return vec;
}

/*---------------------------------------------------------------------------*/
/**
 * @brief Computes the value of the nodal shape function at a given Gauss point
 * (for a given node within in a given element)
 * This method is generic (coming from PASSMO)
 */
/*---------------------------------------------------------------------------*/
Real ArcaneFemFunctions::CellFEMDispatcher::
getShapeFuncVal(Int16 item_type, Integer inod, Real3 coord)
{
 auto f = m_shapefunc[item_type];
 if (f!=nullptr)
   return f(inod,coord);
 return 0.;
}

/*---------------------------------------------------------------------------*/
/**
* @brief Computes the value of the nodal shape function derivatives
 * at a given Gauss point, along 1, 2 and/or 3 directions depending on space dimension
* (for a given node within in a given element)
* This method is generic (coming from PASSMO)
 */
/*---------------------------------------------------------------------------*/
Real3 ArcaneFemFunctions::CellFEMDispatcher::
getShapeFuncDeriv(Int16 item_type, Integer inod, Real3 ref_coord){
 auto f = m_shapefuncderiv[item_type];
 if (f!=nullptr)
   return f(inod,ref_coord);
 return {};
}

