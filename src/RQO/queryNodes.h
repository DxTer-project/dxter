/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2015, The University of Texas and Bryan Marker

    DxTer is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DxTer is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.               

    You should have received a copy of the GNU General Public License
    along with DxTer.  If not, see <http://www.gnu.org/licenses/>.
*/



#pragma once

#include "layers.h"
#include "rqoBasis.h"

#if DORQO

namespace queryNodes
{
    class AndNode;
    class ClauseNode;

    class OrNode
    {
        vector<queryNodes::AndNode> children;

      public:
        OrNode() {};
        virtual void addAnd(AndNode* child) {children.push_back(*child);}
        virtual AndNode* deleteAnd(AndNode* child);
        virtual bool evaluate();
    };

    class AndNode : public OrNode
    {
        vector<queryNodes::ClauseNode> children;

      public:
        AndNode() {};
        virtual void addClause(ClauseNode* child) {children.push_back(*child);}
        virtual ClauseNode* deleteClause(ClauseNode* child);
        virtual bool evaluate();
    };

    class ClauseNode : public AndNode
    {
    public:
        ClauseNode() {};
        virtual bool evaluate();
    };
}





#endif