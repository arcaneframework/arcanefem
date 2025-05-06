//
// Created by ef240508 on 05/05/25.
//

#ifndef ARCANEFEM_MESHTENSORVARIABLEH_H
#define ARCANEFEM_MESHTENSORVARIABLEH_H

#include <arcane/ItemGroup.h>
#include <arcane/VariableScalar.h>
#include <arcane/VariableDataTypeTraits.h>
#include <arcane/VariableTypeInfo.h>
#include <arcane/VariableBuildInfo.h>
#include <arcane/VariableInfo.h>
#include <arcane/VariableFactoryRegisterer.h>

#include "FemUtils.h"
#include "MeshTensorVariable.h"

namespace Arcane
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
template<class ItemType> VariableFactoryRegisterer
MeshVariableScalarTensor2<ItemType>::
m_auto_registerer(_autoCreate,_buildVariableTypeInfo());

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename ItemType> VariableTypeInfo
MeshVariableScalarTensor2<ItemType>::
_buildVariableTypeInfo()
{
  eItemKind ik = ItemTraitsT<ItemType>::kind();
  return BaseClass::_buildVariableTypeInfo(ik);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename ItemType> VariableInfo
MeshVariableScalarTensor2<ItemType>::_buildVariableInfo(const VariableBuildInfo& vbi)
{
  eItemKind ik = ItemTraitsT<ItemType>::kind();
  return BaseClass::_buildVariableInfo(vbi,ik);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType> VariableRef*
MeshVariableScalarTensor2<ItemType>::
_autoCreate(const VariableBuildInfo& vb)
{
  return new ThatClass(vb);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableScalarTensor2<ItemType>::
MeshVariableScalarTensor2(IVariable* var)
: MeshVariableScalarRefT<ItemType, DataType>(var)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableScalarTensor2<ItemType>::
MeshVariableScalarTensor2(const VariableBuildInfo& vb)
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
MeshVariableScalarTensor2<ItemType>::
MeshVariableScalarTensor2(const MeshVariableScalarTensor2<ItemType>& rhs)
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
MeshVariableScalarTensor2<ItemType>::
refersTo(const MeshVariableScalarTensor2<ItemType>& rhs)
{
  MeshVariableScalarRefT<ItemType,DataType>::operator=(rhs);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
typename Arcane::MeshVariableScalarTensor2<ItemType>::GroupType
MeshVariableScalarTensor2<ItemType>::itemGroup() const
{
  return GroupType(this->m_private_part->itemGroup());
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType> VariableFactoryRegisterer
MeshVariableArrayTensor2<ItemType>::
m_auto_registerer(_autoCreate,_buildVariableTypeInfo());

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename ItemType> VariableTypeInfo
MeshVariableArrayTensor2<ItemType>::
_buildVariableTypeInfo()
{
  eItemKind ik = ItemTraitsT<ItemType>::kind();
  eDataType dt = VariableDataTypeTraitsT<DataType>::type();
  return {ik,dt,2,0,false};
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename ItemType> VariableInfo
MeshVariableArrayTensor2<ItemType>::
_buildVariableInfo(const VariableBuildInfo& vbi)
{
  VariableTypeInfo vti = _buildVariableTypeInfo();
  DataStorageTypeInfo sti = vti._internalDefaultDataStorage();
  return {vbi.name(),vbi.itemFamilyName(),vbi.itemGroupName(),vbi.meshName(),vti,sti};
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType> VariableRef*
MeshVariableArrayTensor2<ItemType>::
_autoCreate(const VariableBuildInfo& vb)
{
  return new ThatClass(vb);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<class ItemType>
MeshVariableArrayTensor2<ItemType>::
MeshVariableArrayTensor2(const VariableBuildInfo& vb)
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
MeshVariableArrayTensor2<ItemType>::
MeshVariableArrayTensor2(IVariable* var)
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
MeshVariableArrayTensor2<ItemType>::
MeshVariableArrayTensor2(const MeshVariableArrayTensor2<ItemType>& rhs)
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
MeshVariableArrayTensor2<ItemType>::
refersTo(const MeshVariableArrayTensor2<ItemType>& rhs)
{
  MeshVariableArrayRefT<ItemType, DataType>::operator=(rhs);
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Échange les valeurs de la variable \a rhs avec celles de l'instance.
 *
 * Cette méthode est optimisée pour éviter les recopie et donc l'échange
 * se fait en temps constant. Les variables échangées doivent avoir le
 * même maillage, la même famille et le même groupe. Elles doivent aussi
 * être allouées (IVariable::setUsed()==true)
 */
template<class ItemType,class DataTypeT> void
MeshVariableArrayRefT<ItemType,DataTypeT>::
swapValues(MeshVariableArrayRefT<ItemType,DataType>& rhs)
{
  this->m_private_part->swapValues(*(rhs.m_private_part));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
} // namespace Arcane

#endif //ARCANEFEM_MESHTENSORVARIABLEH_H
