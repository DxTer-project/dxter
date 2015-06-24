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


#include "rqoSort.h"


#if DORQO

Sort::Sort(string sortBy)
  : m_sortBy(sortBy)
{
    static int num = 1;
    m_name = "sort" + std::to_string(num);
    ++num;
}

NodeType Sort::GetType() const
{
  return m_sortBy;
}

void Sort::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const Sort *proj = (Sort*)orig;
  m_name = proj->m_name;
  m_sortBy = proj->m_sortBy;
  Node::Duplicate(orig, shallow, possMerging);
}

const DataTypeInfo& Sort::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void Sort::ClearDataTypeCache()
{
  
}

void Sort::BuildDataTypeCache()
{
  m_dataTypeInfo = InputDataType(0);
  m_dataTypeInfo.m_sortedBy = m_sortBy;
}

void Sort::Prop()
{
  if (m_inputs.size() != 1)
    throw;
  if (m_name.empty())
    throw;

  const DataTypeInfo &in = InputDataType(0);

  if(m_sortBy.empty()) {
    cout << "No sort field for sort\n";
    throw;
  }  

  if (in.m_fields.find(m_sortBy) == in.m_fields.end())
    {
      cout << "sort by is not in input\n";
      throw;
    }
}

void Sort::PrintCode(IndStream &out)
{
  out.Indent();
  string in = GetInputNameStr(0);
  *out << m_name << " = Sort( " << m_sortBy << ", " << in;
  *out << " );\n";
}

Name Sort::GetName(ConnNum num) const
{
  if (num != 0)
    throw;
  Name name(m_name);
  return name;
}

bool RemoveExtraSort::CanApply(const Node *node) const
{
  if (node->m_inputs.size() != 1)
    throw;
  const Node *input = node->Input(0);
  if (input->GetNodeClass() == Sort::GetClass()) {
    return true;
  }
  else
    return false;
}

void RemoveExtraSort::Apply(Node *node) const
{
  Node *input = node->Input(0);
  node->ChangeInput2Way(input, node->InputConnNum(0),
            input->Input(0), input->InputConnNum(0));
  if (input->m_children.empty()) {
    input->m_poss->DeleteChildAndCleanUp(input);
  }
}

bool RemoveRedundantSortBy::CanApply(const Node *node) const
{
  if (node->m_inputs.size() != 1)
    throw;
  const Node *input = node->Input(0);
  if (!node->m_sortBy.empty() && !input->m_sortBy.empty()) {
    return true;
  }
  else
    return false;
}

void RemoveRedundantSortBy::Apply(Node *node) const
{
  Node *input = node->Input(0);
  input->m_sortBy == "";
}

#endif
