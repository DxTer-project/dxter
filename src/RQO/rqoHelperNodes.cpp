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



#include "rqoHelperNodes.h"

#if DORQO


InputNode::InputNode()
:
m_type("InputNode")
{
}


InputNode::InputNode(string name, string sortBy, set<string> fields, string fileName, string query)
:
  m_type(name),
  m_varName(name),
  m_fileName(fileName),
  m_query(query)
{
  m_dataTypeInfo.m_sortedBy = sortBy;
  m_dataTypeInfo.m_fields = fields;
}

void InputNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  Node::Duplicate(orig, shallow, possMerging);
  const InputNode *node = (InputNode*)orig;
  m_type = node->m_type;
  m_dataTypeInfo = node->m_dataTypeInfo;
  m_varName = node->m_varName;
  m_fileName = node->m_fileName;
  m_query = node->m_query;
}

const DataTypeInfo& InputNode::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void InputNode::Prop()
{
  if (!m_inputs.empty())
    throw;
  if (m_type.empty() || m_varName.empty())
    throw;
  if (m_dataTypeInfo.m_fields.find(m_dataTypeInfo.m_sortedBy)
      == m_dataTypeInfo.m_fields.end())
    throw;
  if(m_fileName.empty() || m_query.empty())
    throw;
  //add something to check for valid queries or valid files
}

void InputNode::PrintCode(IndStream &out)
{
  out.Indent();
  *out << m_varName << " = Scan( " << m_fileName << ", " << m_query << ");\n";
}

Name InputNode::GetName(ConnNum num) const
{
  if (num != 0)
    throw;
  return m_varName;
}

const DataTypeInfo& OutputNode::DataType(ConnNum num) const
{
  throw;
}

void OutputNode::Prop()
{
  if (m_inputs.empty())
    throw;
  if (!m_children.empty())
    throw;
}

void OutputNode::PrintCode(IndStream &out)
{

}
#endif //DORQO

