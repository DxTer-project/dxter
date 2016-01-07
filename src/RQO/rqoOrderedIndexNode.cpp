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



#include "rqoOrderedIndexNode.h"

#if DORQO


OrderedIndexedNode::OrderedIndexedNode()
:
  InputNode()
{
}


OrderedIndexedNode::OrderedIndexedNode(string name, string sortBy, set<string> fields, string fileName, string query, int index)
:
  m_index(index)
{
  InputNode(name, sortBy, fields, fileName, query);
  m_varName = name;
  m_query = query;
  m_fileName = fileName;
  m_dataTypeInfo.m_sortedBy = sortBy;
  m_dataTypeInfo.m_fields = fields;
}

void OrderedIndexedNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  Node::Duplicate(orig, shallow, possMerging);
  const OrderedIndexedNode *node = (OrderedIndexedNode*)orig;
  m_type = node->m_type;
  m_dataTypeInfo = node->m_dataTypeInfo;
  m_varName = node->m_varName;
  m_fileName = node->m_fileName;
  m_query = node->m_query;
  m_relation = node->m_relation;
  m_index = node->m_index;
}

const DataTypeInfo& OrderedIndexedNode::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void OrderedIndexedNode::Prop()
{
  InputNode::Prop();
}

void OrderedIndexedNode::PrintCode(IndStream &out)
{
  out.Indent();
  *out << m_varName << " = orderedindexFunc( " << m_fileName << ", " << m_query << ", " << m_index;
  for(auto string : m_dataTypeInfo.m_fields)
  {
    *out << "," << string;
  }

  *out << ");\n";
}

Name OrderedIndexedNode::GetName(ConnNum num) const
{
  if (num != 0)
    throw;
  return m_varName;
}

#endif