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

#include "rqoProj.h"


#if DORQO

Projection::Projection(string sortBy,
    set<string> &inFields)
  : m_inFields(inFields),
    m_sortBy(sortBy)
{
    static int num = 1;
    m_name = "projection" + std::to_string(num);
    ++num;
}

NodeType Projection::GetType() const
{
  string ret = m_sortBy;
  set<string>::const_iterator iter = m_inFields.begin();
  for(; iter != m_inFields.end(); ++iter) {
    ret += "," + *iter;
  }
  return ret;
}

void Projection::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const Projection *proj = (Projection*)orig;
  m_name = proj->m_name;
  m_sortBy = proj->m_sortBy;
  m_inFields = proj->m_inFields;
  Node::Duplicate(orig, shallow, possMerging);
}

const DataTypeInfo& Projection::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void Projection::ClearDataTypeCache()
{
  
}

void Projection::BuildDataTypeCache()
{
  m_dataTypeInfo.m_sortedBy = m_sortBy;
  m_dataTypeInfo.m_fields = m_inFields;
}

void Projection::Prop()
{
  if (m_inputs.size() != 1)
    throw;
  if (m_inFields.empty())
    throw;
  if (m_sortBy.empty())
    throw;
  if (m_name.empty())
    throw;

  const DataTypeInfo &in = InputDataType(0);


  for (auto str : m_inFields) {
    if (in.m_fields.find(str) == in.m_fields.end())
      throw;
  }

  if (m_inFields.find(m_sortBy) == m_inFields.end())
    {
      cout << "sort by is not in input\n";
      throw;
    }
}

void Projection::PrintCode(IndStream &out)
{
  out.Indent();
  string in = GetInputNameStr(0);
  *out << m_name << " = Projection( " << m_sortBy << ", " << in;
  set<string>::iterator iter = m_inFields.begin();
  for(; iter != m_inFields.end(); ++iter) {
    *out << ", " << in << "." << *iter;
  }
  *out << " );\n";
}

Name Projection::GetName(ConnNum num) const
{
  if (num != 0)
    throw;
  Name name(m_name);
  return name;
}


bool RemoveExtraProjection::CanApply(const Node *node) const
{
  if (node->m_inputs.size() != 1)
    throw;
  const Node *input = node->Input(0);
  if (input->GetNodeClass() == Projection::GetClass()) {
    return true;
  }
  else
    return false;
}

void RemoveExtraProjection::Apply(Node *node) const
{
  Node *input = node->Input(0);
  node->ChangeInput2Way(input, node->InputConnNum(0),
			input->Input(0), input->InputConnNum(0));
  if (input->m_children.empty()) {
    input->m_poss->DeleteChildAndCleanUp(input);
  }
}

#endif
