//
// Created by evely on 03/12/2025.
//

#ifndef GAUSSITEM_H
#define GAUSSITEM_H

#include "arcane/Item.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class ARCANE_CORE_EXPORT GaussPoint
: public Item
{
  using ThatClass = GaussPoint;
  // Pour accéder aux constructeurs privés
  friend class ItemEnumeratorBaseT<ThatClass>;
  friend class ItemConnectedEnumeratorBaseT<ThatClass>;
  friend class ItemVectorT<ThatClass>;
  friend class ItemVectorViewT<ThatClass>;
  friend class ItemConnectedListViewT<ThatClass>;
  friend class ItemVectorViewConstIteratorT<ThatClass>;
  friend class ItemConnectedListViewConstIteratorT<ThatClass>;
  friend class SimdItemT<ThatClass>;
  friend class ItemInfoListViewT<ThatClass>;
  friend class ItemLocalIdToItemConverterT<ThatClass>;

 public:

  class ARCANE_DEPRECATED_REASON("Y2024: Use GaussLocalId instead") Index
  : public Item::Index
  {
   public:
    typedef Item::Index Base;
   public:
    explicit Index(Int32 id) : Base(id){}
    Index(GaussPoint item) : Base(item){}
    operator GaussLocalId() const { return GaussLocalId{localId()}; }
  };

 protected:
  constexpr GaussPoint(Int32 local_id,ItemSharedInfo* shared_info)
  : Item(local_id,shared_info) {}

 public:
  typedef GaussLocalId LocalIdType;
  GaussPoint() = default;
  GaussPoint(ItemInternal* ainternal) : Item(ainternal)
  { ARCANE_CHECK_KIND(isGauss); }
  constexpr GaussPoint(const ItemBase& abase) : Item(abase)
  { ARCANE_CHECK_KIND(isGauss); }
  constexpr explicit GaussPoint(const Item& aitem) : Item(aitem)
  { ARCANE_CHECK_KIND(isGauss); }
  GaussPoint(const ItemInternalPtr* internals,Int32 local_id) : Item(internals,local_id)
  { ARCANE_CHECK_KIND(isGauss); }
  GaussPoint& operator=(ItemInternal* ainternal)
  {
    _set(ainternal);
    return (*this);
  }

 public:
  constexpr eItemKind kind() const { return IK_Gauss; }
  constexpr GaussLocalId itemLocalId() const { return GaussLocalId{ m_local_id }; }
  constexpr Int32 nbEdge() const { return _nbEdge(); }
  constexpr Int32 nbFace() const { return _nbFace(); }
  Int32 nbCell() const { return _nbCell(); }
  inline Edge edge(Int32 i) const;
  inline Face face(Int32 i) const;
  inline Cell cell(Int32 i) const;
  EdgeLocalId edgeId(Int32 i) const { return _edgeId(i); }
  FaceLocalId faceId(Int32 i) const { return _faceId(i); }
  CellLocalId cellId(Int32 i) const { return _cellId(i); }
  EdgeConnectedListViewType edges() const { return _edgeList(); }
  FaceConnectedListViewType faces() const { return _faceList(); }
  CellConnectedListViewType cells() const { return _cellList(); }
  EdgeLocalIdView edgeIds() const { return _edgeIds(); }
  FaceLocalIdView faceIds() const { return _faceIds(); }
  CellLocalIdView cellIds() const { return _cellIds(); }

  // AMR
  CellVectorView _internalActiveCells(Int32Array& local_ids) const
  {
    return _toItemBase()._internalActiveCells2(local_ids);
  }

  ARCANE_DEPRECATED_REASON("Y2022: Do not use this operator. Use operator '.' instead")
  GaussPoint* operator->() { return this; }

  ARCANE_DEPRECATED_REASON("Y2022: Do not use this operator. Use operator '.' instead")
  const GaussPoint* operator->() const { return this; }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

constexpr inline GaussPoint Item::
_gausspoint(Int32 index) const
{
  return GaussPoint(_connectivity()->nodeBase(m_local_id,index));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
class ARCANE_CORE_EXPORT ItemWithNodes
: public Item
{
  using ThatClass = ItemWithNodes;
  // Pour accéder aux constructeurs privés
  friend class ItemEnumeratorBaseT<ThatClass>;
  friend class ItemConnectedEnumeratorBaseT<ThatClass>;
  friend class ItemVectorT<ThatClass>;
  friend class ItemVectorViewT<ThatClass>;
  friend class ItemConnectedListViewT<ThatClass>;
  friend class ItemVectorViewConstIteratorT<ThatClass>;
  friend class ItemConnectedListViewConstIteratorT<ThatClass>;
  friend class SimdItemT<ThatClass>;
  friend class ItemInfoListViewT<ThatClass>;
  friend class ItemLocalIdToItemConverterT<ThatClass>;

 protected:
  constexpr ItemWithNodes(Int32 local_id,ItemSharedInfo* shared_info)
  : Item(local_id,shared_info) {}

 public:

  ItemWithNodes() = default;
  ItemWithNodes(ItemInternal* ainternal) : Item(ainternal)
  { ARCANE_CHECK_KIND(isItemWithNodes); }
  constexpr ItemWithNodes(const ItemBase& abase) : Item(abase)
  { ARCANE_CHECK_KIND(isItemWithNodes); }
  constexpr explicit ItemWithNodes(const Item& aitem) : Item(aitem)
  { ARCANE_CHECK_KIND(isItemWithNodes); }
  ItemWithNodes(const ItemInternalPtr* internals,Int32 local_id)
  : Item(internals,local_id)
  { ARCANE_CHECK_KIND(isItemWithNodes); }
  ItemWithNodes& operator=(ItemInternal* ainternal)
  {
    _set(ainternal);
    return (*this);
  }

 public:
  Int32 nbNode() const { return _nbNode(); }
  Node node(Int32 i) const { return _node(i); }
  NodeConnectedListViewType nodes() const { return _nodeList(); }
  NodeLocalIdView nodeIds() const { return _nodeIds(); }
  NodeLocalId nodeId(Int32 index) const { return _nodeId(index); }
  Int32 nbLinearNode() const { return _nbLinearNode(); }

 public:

  ARCANE_DEPRECATED_REASON("Y2022: Do not use this operator. Use operator '.' instead")
  ItemWithNodes* operator->() { return this; }

  ARCANE_DEPRECATED_REASON("Y2022: Do not use this operator. Use operator '.' instead")
  const ItemWithNodes* operator->() const { return this; }
};

#endif //GAUSSITEM_H
