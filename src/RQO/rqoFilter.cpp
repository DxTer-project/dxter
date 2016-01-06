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

#include "rqoFilter.h"


#if DORQO

Filter::Filter()
  : Sortable()
  {
    
  }

Filter::Filter(string sortBy,
    set<string> &inFields)
  : Sortable(sortBy),
    m_inFields(inFields)
{
    static int num = 1;
    m_name = "filter" + std::to_string(num);
    ++num;
}

NodeType Filter::GetType() const
{
  string ret = m_sortBy;
  set<string>::const_iterator iter = m_inFields.begin();
  for(; iter != m_inFields.end(); ++iter) {
    ret += "," + *iter;
  }
  return ret;
}

void Filter::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const Filter *fil = (Filter*)orig;
  m_name = fil->m_name;
  m_sortBy = fil->m_sortBy;
  m_inFields = fil->m_inFields;
  Node::Duplicate(orig, shallow, possMerging);
}

const DataTypeInfo& Filter::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void Filter::ClearDataTypeCache()
{
  
}

void Filter::BuildDataTypeCache()
{
  m_dataTypeInfo.m_sortedBy = m_sortBy;
  m_dataTypeInfo.m_fields = m_inFields;
}

void Filter::Prop()
{
  if (m_inputs.size() != 1)
    throw;
  if (m_inFields.empty())
    throw;
  if (m_name.empty())
    throw;

  const DataTypeInfo &in = InputDataType(0);


  for (auto str : m_inFields) {
    if (in.m_fields.find(str) == in.m_fields.end())
      throw;
  }


}

void Filter::PrintCode(IndStream &out)
{
  out.Indent();
  string in = GetInputNameStr(0);
  *out << m_name << " = filter(" << m_sortBy << "," << in;
  set<string>::iterator iter = m_inFields.begin();
  for(; iter != m_inFields.end(); ++iter) {
    *out << "," << in << "." << *iter;
  }
  *out << ");\n";
}


#endif
