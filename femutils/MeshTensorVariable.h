//
// Created by ef240508 on 05/05/25.
//

#ifndef ARCANEFEM_MESHTENSORVARIABLE_H
#define ARCANEFEM_MESHTENSORVARIABLE_H

#include <arcane/MeshVariableRef.h>
#include <arcane/MeshVariableInfo.h>
#include <arcane/PrivateVariableScalar.h>
#include <arcane/ItemEnumerator.h>
#include <arcane/ItemGroupRangeIterator.h>
#include <arcane/ItemPairEnumerator.h>
#include <arcane/core/GroupIndexTable.h>

#include "FemUtils.h"

namespace Arcane
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
template<class ItemTypeT>
class MeshVariableScalarTensor2:
public MeshVariableScalarRefT<ItemTypeT,FemUtils::Tensor2>
{
  typedef FemUtils::Tensor2 DataType;
  typedef ItemTypeT ItemType;
  typedef UniqueArray<DataType> ValueType;
  typedef const DataType& ConstReturnReferenceType;
  typedef DataType& ReturnReferenceType;

 protected:

  typedef MeshVariableScalarRefT<ItemTypeT,DataType> BaseClass;

  typedef typename ItemType::Index ItemIndexType;
  typedef typename ItemType::LocalIdType ItemLocalIdType;
  typedef typename ItemTraitsT<ItemType>::ItemGroupType GroupType;
  typedef MeshVariableScalarTensor2<ItemType> ThatClass;
  typedef typename BaseClass::DataTypeReturnReference DataTypeReturnReference;

 public:

  //! Construit une référence à la variable spécifiée dans \a vb
  explicit ARCANE_CORE_EXPORT MeshVariableScalarTensor2(const VariableBuildInfo& vb);

  //! Construit une référence à partir de \a var
  explicit ARCANE_CORE_EXPORT MeshVariableScalarTensor2(IVariable* var);

  //! Construit une référence à partir de \a rhs
  ARCANE_CORE_EXPORT MeshVariableScalarTensor2(const MeshVariableScalarTensor2<ItemType>& rhs);

  //! Positionne la référence de l'instance à la variable \a rhs.
  ARCANE_CORE_EXPORT void refersTo(const MeshVariableScalarTensor2<ItemType>& rhs);

  void resize(Int32 dim);

//  ThatClass& operator=(const ThatClass& rhs) = delete;
  //! Constructeur vide
  MeshVariableScalarTensor2()= default;

  //! Groupe associé à la grandeur
  ARCANE_CORE_EXPORT GroupType itemGroup() const;

  //! Valeur non modifiable de l'entité \a item
  const DataType& operator[](ItemLocalIdType i) const { return this->_value(i.localId()); }

  //! Valeur modifiable de l'entité \a item
  DataTypeReturnReference operator[](ItemLocalIdType i) { return this->_value(i.localId()); }

  //! Valeur non modifiable de l'entité \a item
  const DataType& operator()(ItemLocalIdType i) const { return this->_value(i.localId()); }

  //! Valeur modifiable de l'entité \a item
  DataTypeReturnReference operator()(ItemLocalIdType i) { return this->_value(i.localId()); }

  //! Valeur modifiable de l'entité \a item
  DataType& getReference(ItemLocalIdType item)
  {
    return this->_value(item.localId());
  }

  const DataType& item(ItemLocalIdType i) const
  {
    return this->_value(i.localId());
  }
  void setItem(ItemLocalIdType i,const DataType& v)
  {
    this->_value(i.localId()) = v;
  }

 private:

  static VariableFactoryRegisterer m_auto_registerer;
  static VariableInfo _buildVariableInfo(const VariableBuildInfo& vbi);
  static VariableTypeInfo _buildVariableTypeInfo();
  static VariableRef* _autoCreate(const VariableBuildInfo& vb);

};

template<class ItemTypeT>
class MeshVariableArrayTensor2:
public MeshVariableArrayRefT<ItemTypeT,FemUtils::Tensor2>
{
 public:

  typedef FemUtils::Tensor2 DataType;
  typedef ItemTypeT ItemType;
  typedef UniqueArray2<DataType> ValueType;
  typedef ConstArrayView<DataType> ConstReturnReferenceType;
  typedef ArrayView<DataType> ReturnReferenceType;

 protected:

  typedef MeshVariableArrayRefT<ItemTypeT,DataType> BaseClass;
  typedef typename ItemTraitsT<ItemType>::ItemGroupType GroupType;
  typedef MeshVariableArrayRefT<ItemType,DataType> ThatClass;

  typedef Array2View<DataType> ArrayBase;
  typedef ArrayView<DataType> ArrayType;
  typedef ConstArrayView<DataType> ConstArrayType;
  typedef Array2VariableT<DataType> PrivatePartType;

 public:
  //! Constructeur vide
  MeshVariableArrayTensor2()= default;

  //! Construit une référence à la variable spécifiée dans \a vb
  explicit ARCANE_CORE_EXPORT MeshVariableArrayTensor2(const VariableBuildInfo& b);
  //! Construit une référence à partir de \a var
  explicit ARCANE_CORE_EXPORT MeshVariableArrayTensor2(IVariable* var);
  //! Construit une référence à partir de \a rhs
  ARCANE_CORE_EXPORT MeshVariableArrayTensor2(const MeshVariableArrayTensor2<ItemType>& rhs);
  //! Positionne la référence de l'instance à la variable \a rhs.
  ARCANE_CORE_EXPORT void refersTo(const MeshVariableArrayTensor2<ItemType>& rhs);

 public:

  private:

  static VariableFactoryRegisterer m_auto_registerer;
  static VariableInfo _buildVariableInfo(const VariableBuildInfo& vbi);
  static VariableTypeInfo _buildVariableTypeInfo();
  static VariableRef* _autoCreate(const VariableBuildInfo& vb);
};

/*!
 * \ingroup Variable
 * \brief Classe de base d'une variable type Tensor2 sur des entités du maillage.
 */
template<>
class MeshVariableInfoT<Node,Arcane::FemUtils::Tensor2,0>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableScalarRefT<Node,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef VariableArrayT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<Node,Arcane::FemUtils::Tensor2,1>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableArrayRefT<Node,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef Array2VariableT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<Edge,Arcane::FemUtils::Tensor2,0>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableScalarRefT<Edge,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef VariableArrayT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<Edge,Arcane::FemUtils::Tensor2,1>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableArrayRefT<Edge,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef Array2VariableT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<Face,Arcane::FemUtils::Tensor2,0>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableScalarRefT<Node,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef VariableArrayT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<Face,Arcane::FemUtils::Tensor2,1>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableArrayRefT<Node,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef Array2VariableT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<Cell,Arcane::FemUtils::Tensor2,0>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableScalarRefT<Cell,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef VariableArrayT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<Cell,Arcane::FemUtils::Tensor2,1>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableArrayRefT<Cell,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef Array2VariableT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<DoF,Arcane::FemUtils::Tensor2,0>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableScalarRefT<DoF,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef VariableArrayT< Arcane::FemUtils::Tensor2 > PrivateType;
};

template<>
class MeshVariableInfoT<DoF,Arcane::FemUtils::Tensor2,1>
{
 public:
  //! Type de la référence  la variable
  typedef MeshVariableArrayRefT<DoF,Arcane::FemUtils::Tensor2> RefType;
  //! Type de la partie privé de la variable
  typedef Array2VariableT< Arcane::FemUtils::Tensor2 > PrivateType;
};

typedef MeshVariableArrayTensor2<Arcane::DoF> VariableDoFArrayTensor2;

} // namespace Arcane

#endif //ARCANEFEM_MESHTENSORVARIABLE_H
