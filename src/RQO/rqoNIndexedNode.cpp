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



#include "rqoNIndexedNode.h"

#if DORQO


NIndexedNode::NIndexedNode()
:
  InputNode()
{
}


NIndexedNode::NIndexedNode(string name, string sortBy, set<string> fields, string fileName, string query, set<int> indeces)
:
  m_indeces(indeces)
{
  InputNode(name, sortBy, fields, fileName, query);
  m_varName = name;
  m_query = query;
  m_fileName = fileName;
  m_dataTypeInfo.m_sortedBy = sortBy;
  m_dataTypeInfo.m_fields = fields;
}

void NIndexedNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  Node::Duplicate(orig, shallow, possMerging);
  const NIndexedNode *node = (NIndexedNode*)orig;
  m_type = node->m_type;
  m_dataTypeInfo = node->m_dataTypeInfo;
  m_varName = node->m_varName;
  m_fileName = node->m_fileName;
  m_query = node->m_query;
  m_relation = node->m_relation;
  m_indeces = node->m_indeces;
}

const DataTypeInfo& NIndexedNode::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void NIndexedNode::Prop()
{
  InputNode::Prop();
}

void NIndexedNode::PrintCode(IndStream &out)
{
  out.Indent();
  *out << m_varName << " = nindexFunc(" << m_fileName << "," << m_query << ",[";

  set<int>::iterator iter = m_indeces.begin();
  for(; iter != m_indeces.end(); ++iter)
  {
    *out << (*iter) << ",";
  }


  *out << "]);\n";
}

Name NIndexedNode::GetName(ConnNum num) const
{
  if (num != 0)
    throw;
  return m_varName;
}

#endif