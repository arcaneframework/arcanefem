//
// Created by aa286863 on 16/06/2026.
//

#ifndef ARCANEFEM_MESHMATRIXVARIABLESHELPERS_H
#define ARCANEFEM_MESHMATRIXVARIABLESHELPERS_H

#include <arcane/ItemGroup.h>
#include <arcane/VariableScalar.h>
#include <arcane/VariableDataTypeTraits.h>
#include <arcane/VariableTypeInfo.h>
#include <arcane/VariableBuildInfo.h>
#include <arcane/VariableInfo.h>
#include <arcane/VariableFactoryRegisterer.h>

#include "FemUtils.h"
#include "MeshMatrixVariables.h"

namespace Arcane
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
template<class ItemType> VariableFactoryRegisterer
MeshVariableScalarRealMatrix3x3<ItemType>::
m_auto_registerer(_autoCreate,_buildVariableTypeInfo());

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename ItemType> VariableTypeInfo
MeshVariableScalarRealMatrix3x3<ItemType>::
_buildVariableTypeInfo()
{
  eItemKind ik = ItemTraitsT<ItemType>::kind();
  return BaseClass::_buildVariableTypeInfo(ik);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename ItemType> VariableInfo
MeshVariableScalarRealMatrix3x3<ItemType>::_buildVariableInfo(const VariableBuildInfo& vbi)
{
  eItemKind ik = ItemTraitsT<ItemType>::kind();
  return BaseClass::_buildVariableInfo(vbi,ik);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType> VariableRef*
MeshVariableScalarRealMatrix3x3<ItemType>::
_autoCreate(const VariableBuildInfo& vb)
{
  return new ThatClass(vb);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableScalarRealMatrix3x3<ItemType>::
MeshVariableScalarRealMatrix3x3(IVariable* var)
: MeshVariableScalarRefT<ItemType, DataType>(var)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableScalarRealMatrix3x3<ItemType>::
MeshVariableScalarRealMatrix3x3(const VariableBuildInfo& vb)
: MeshVariableScalarRefT<ItemType, DataType>(vb,ItemTraitsT<ItemType>::kind())
{
  // Normalement, c'est à cette classe de faire l'initilisation mais
  // comme cette classe est juste un wrapper autour de ItemVariableArrayRefT
  // et ne fait rien d'autre, on laisse l'initialisation à la classe de base,
  // ce qui permet de fabriquer de manière générique une variable sur
  // une entité du maillage à partir de son genre.
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableScalarRealMatrix3x3<ItemType>::
MeshVariableScalarRealMatrix3x3(const MeshVariableScalarRealMatrix3x3<ItemType>& rhs)
: MeshVariableScalarRefT<ItemType,DataType>(rhs)
{
  // Normalement, c'est à cette classe de faire l'initilisation mais
  // comme cette classe est juste un wrapper autour de ItemVariableArrayRefT
  // et ne fait rien d'autre, on laisse l'initialisation à la classe de base,
  // ce qui permet de fabriquer de manière générique une variable sur
  // une entité du maillage à partir de son genre.
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType> void
MeshVariableScalarRealMatrix3x3<ItemType>::
refersTo(const MeshVariableScalarRealMatrix3x3<ItemType>& rhs)
{
  MeshVariableScalarRefT<ItemType,DataType>::operator=(rhs);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
typename Arcane::MeshVariableScalarRealMatrix3x3<ItemType>::GroupType
MeshVariableScalarRealMatrix3x3<ItemType>::itemGroup() const
{
  return GroupType(this->m_private_part->itemGroup());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType> VariableFactoryRegisterer
MeshVariableArrayRealMatrix3x3<ItemType>::
m_auto_registerer(_autoCreate,_buildVariableTypeInfo());

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename ItemType> VariableTypeInfo
MeshVariableArrayRealMatrix3x3<ItemType>::
_buildVariableTypeInfo()
{
  eItemKind ik = ItemTraitsT<ItemType>::kind();
  eDataType dt = VariableDataTypeTraitsT<DataType>::type();
  return {ik,dt,2,0,false};
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename ItemType> VariableInfo
MeshVariableArrayRealMatrix3x3<ItemType>::
_buildVariableInfo(const VariableBuildInfo& vbi)
{
  VariableTypeInfo vti = _buildVariableTypeInfo();
  DataStorageTypeInfo sti = vti._internalDefaultDataStorage();
  return {vbi.name(),vbi.itemFamilyName(),vbi.itemGroupName(),vbi.meshName(),vti,sti};
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType> VariableRef*
MeshVariableArrayRealMatrix3x3<ItemType>::
_autoCreate(const VariableBuildInfo& vb)
{
  return new ThatClass(vb);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableArrayRealMatrix3x3<ItemType>::
MeshVariableArrayRealMatrix3x3(const VariableBuildInfo& vb)
: MeshVariableArrayRefT<ItemType, DataType>(vb,ItemTraitsT<ItemType>::kind())
{
  // Normalement, c'est à cette classe de faire l'initilisation mais
  // comme cette classe est juste un wrapper autour de ItemVariableArrayRefT
  // et ne fait rien d'autre, on laisse l'initialisation à la classe de base,
  // ce qui permet de fabriquer de manière générique une variable sur
  // une entité du maillage à partir de son genre.
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableArrayRealMatrix3x3<ItemType>::
MeshVariableArrayRealMatrix3x3(IVariable* var)
: MeshVariableArrayRefT<ItemType, DataType>(var)
{
  // Normalement, c'est à cette classe de faire l'initilisation mais
  // comme cette classe est juste un wrapper autour de ItemVariableArrayRefT
  // et ne fait rien d'autre, on laisse l'initialisation à la classe de base,
  // ce qui permet de fabriquer de manière générique une variable sur
  // une entité du maillage à partir de son genre.
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableArrayRealMatrix3x3<ItemType>::
MeshVariableArrayRealMatrix3x3(const MeshVariableArrayRealMatrix3x3<ItemType>& rhs)
: MeshVariableArrayRefT<ItemType, DataType>(rhs)
{
  // Normalement, c'est à cette classe de faire l'initilisation mais
  // comme cette classe est juste un wrapper autour de ItemVariableArrayRefT
  // et ne fait rien d'autre, on laisse l'initialisation à la classe de base,
  // ce qui permet de fabriquer de manière générique une variable sur
  // une entité du maillage à partir de son genre.
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
void
MeshVariableArrayRealMatrix3x3<ItemType>::
refersTo(const MeshVariableArrayRealMatrix3x3<ItemType>& rhs)
{
  MeshVariableArrayRefT<ItemType, DataType>::operator=(rhs);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
} // namespace Arcane

#endif //ARCANEFEM_MESHMATRIXVARIABLESHELPERS_H
