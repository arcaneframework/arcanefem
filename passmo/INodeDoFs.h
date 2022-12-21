//
// Created by EF on 15/12/22.
//

#ifndef PASSMO_INODEDOFS_H
#define PASSMO_INODEDOFS_H

#include <arcane/ItemTypes.h>

using namespace Arcane;

/**
 * Interface  of NodeDoFs service for FEM modules
 */
class INodeDoFs
{
public:
    INodeDoFs() {};
    virtual ~INodeDoFs() {};

public:
    virtual Integer addDoF(const Integer& begin_uid) = 0;
    virtual Integer removeDoF(const Integer& begin_lid) = 0;
    virtual Integer addDoFs(const Integer& begin_uid, const Integer& size) = 0;
    virtual Integer removeDoFs(const Integer& begin_lid, const Integer& size) = 0;
    virtual DoFGroup createDoFGroup(const Integer& begin_lid, const Integer& size) = 0;
    virtual VariableDoFReal doFVariable() = 0;
    virtual VariableDoFArrayReal doFVectorVariable(const Integer& size) = 0;
};

#endif //PASSMO_INODEDOFS_H
