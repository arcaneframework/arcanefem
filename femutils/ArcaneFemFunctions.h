#ifndef ARCANE_FEM_FUNCTIONS_H
#define ARCANE_FEM_FUNCTIONS_H

#include <arcane/VariableTypes.h>

using namespace Arcane;

class ArcaneFemFunctions
{
public:
  static Real computeAreaTriangle3(Cell cell, VariableNodeReal3& node_coord);
  static Real computeEdgeLength2(Face face, VariableNodeReal3& node_coord);
  static Real2 computeEdgeNormal2(Face face, VariableNodeReal3& node_coord);
};

#endif // ARCANE_FEM_FUNCTIONS_H
