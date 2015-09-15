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
#include "rqoTuple.h"
#include "rqoPredicate.h"

#if DORQO

//The hierarchy of these nodes enforces correctness in queries. 
namespace queryNodes
{
    class AndNode;
    class ClauseNode;


    class OrNode
    {
        string m_id;
        vector<queryNodes::AndNode> children;

      public:
        OrNode();
        //This adds a child node to the vector children.
        virtual void addAnd(AndNode* child) {children.push_back(*child);}
        //This Deletes a specified child node in children
        virtual void deleteAnd(AndNode* child);
        //Evaluate surveys the query and makes sure it is possible. All of the evaluates
        virtual bool evaluate(Tuple tuple);
        //Return the id of the node
        virtual string getId() {return m_id;}
    };

    class AndNode : public OrNode
    {
        string m_id;
        vector<queryNodes::ClauseNode> children;

      public:
        AndNode();
        //Adds a clause node to children.
        virtual void addClause(ClauseNode* child) {children.push_back(*child);}
        //Deletes a specified node from children.
        virtual void deleteClause(ClauseNode* child);
        virtual bool evaluate(Tuple tuple);
    };

    class ClauseNode : public AndNode
    {
        
    public:
        string m_id;

        ClauseNode();
    };

    class FieldValue : public ClauseNode
    {
    public:
        FieldValue();
        virtual bool evaluate(Tuple tuple);
    };

    class FieldField : public ClauseNode
    {
    public:
        FieldField();
        virtual bool evaluate(Tuple tuple);
    };

    class FieldSet : public ClauseNode
    {
    public:
        FieldSet();
        virtual bool evaluate(Tuple tuple);
    };
}





#endif