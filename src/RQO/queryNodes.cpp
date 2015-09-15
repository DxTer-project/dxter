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


#include "queryNodes.h"


#if DORQO

using namespace queryNodes;

OrNode::OrNode()
{
    static int num = 1;
    m_id = "orNode" + std::to_string(num);
    ++num;
}

bool OrNode::evaluate(Tuple tuple)
{
    bool ret = false;
    vector<AndNode>::iterator iter = children.begin();
    for(; iter != children.end(); iter++)
    {
        ret = ret || (*iter).evaluate(tuple);
    }
    return ret;
}


void OrNode::deleteAnd(AndNode* child)
{
    vector<AndNode>::iterator iter = children.begin();
    int index = 0;
    for(; iter != children.end(); iter++)
    {
        if((*iter).getId() == child->getId())
        {
            children.erase(children.begin() + index);
            break;
        }
        index++;
    }
}

AndNode::AndNode()
{
    static int num = 1;
    m_id = "andNode" + std::to_string(num);
    ++num;
}

void AndNode::deleteClause(ClauseNode* child)
{
    vector<ClauseNode>::iterator iter = children.begin();
    int index = 0;
    for(; iter != children.end(); iter++)
    {
        if((*iter).getId() == child->getId())
        {
            children.erase(children.begin() + index);
            break;
        }
        index++;
    }
}

bool AndNode::evaluate(Tuple tuple)
{
    bool ret = false;
    vector<ClauseNode>::iterator iter = children.begin();
    for(; iter != children.end(); iter++)
    {
        ret = ret || (*iter).evaluate(tuple);
    }
    return ret;
}

ClauseNode::ClauseNode()
{
    static int num = 1;
    m_id = "clauseNode" + std::to_string(num);
    ++num;
}

FieldValue::FieldValue()
{
    static int num = 1;
    m_id = "fieldValue" + std::to_string(num);
    ++num;
}

FieldField::FieldField()
{
    static int num = 1;
    m_id = "fieldField" + std::to_string(num);
    ++num;
}

FieldSet::FieldSet()
{
    static int num = 1;
    m_id = "fieldSet" + std::to_string(num);
    ++num;
}

bool FieldValue::evaluate(Tuple tuple)
{
    bool ret = false;
    return ret;
}

bool FieldField::evaluate(Tuple tuple)
{
    bool ret = false;
    return ret;
}

bool FieldSet::evaluate(Tuple tuple)
{
    bool ret = false;
    return ret;
}

#endif