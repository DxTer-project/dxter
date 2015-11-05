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


#include "rqoUnion.h"


#if DORQO

Union::Union()
  : Sortable()
{

}

Union::Union(string sortBy, 
	   vector<string> inFields)
  : Sortable(sortBy),
    m_fields(inFields)
{
  static int num = 1;
  m_name = "union" + std::to_string(num);
  ++num;
}

NodeType Union::GetType() const
{
  string ret = m_sortBy;

  vector<string>::const_iterator iter = m_fields.begin();
  for(; iter != m_fields.end(); ++iter) {
    ret += "," + *iter;
  }
  return ret;
}

void Union::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const Union *unionNode = (Union*)orig;
  m_name = unionNode->m_name;
  m_sortBy = unionNode->m_sortBy;
  m_fields = unionNode->m_fields;
  Node::Duplicate(orig, shallow, possMerging);
}

const DataTypeInfo& Union::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void Union::ClearDataTypeCache()
{
  
}

void Union::BuildDataTypeCache()
{
  m_dataTypeInfo.m_sortedBy = m_sortBy;
  m_dataTypeInfo.m_fields = InputDataType(0).m_fields;
  m_dataTypeInfo.m_fields.insert(InputDataType(1).m_fields.begin(),
				 InputDataType(1).m_fields.end());
}

void Union::Prop()
{
  //Node::Prop();
  if (m_inputs.size() != 2)
    throw;
  if (m_fields.empty())
    throw;
  if (m_name.empty())
    throw;

  const DataTypeInfo &in0 = InputDataType(0);
  const DataTypeInfo &in1 = InputDataType(1);


  for (auto str : m_fields) {
    if (in0.m_fields.find(str) == in0.m_fields.end() || 
          in1.m_fields.find(str) == in1.m_fields.end()) {
      cout << "input does not include field over which we are Unioning\n";
      throw;
    }
  }


  if(!m_sortBy.empty())
  {
    if (in0.m_fields.find(m_sortBy) == in0.m_fields.end()
      && in1.m_fields.find(m_sortBy) == in1.m_fields.end())
    {
      cout << "sort by is not in either input\n";
      throw;
    }
  }
  
}

void Union::PrintCode(IndStream &out)
{
  out.Indent();
  string in0 = GetInputNameStr(0);
  string in1 = GetInputNameStr(1);
  *out << m_name << " = unionFunc(" << m_sortBy << ","
    << in0 << "," << in1;
  vector<string>::iterator iter = m_fields.begin(); 
  for(; iter != m_fields.end(); ++iter) {
    *out << "," << *iter;
  }
  *out << ");\n";
}


#endif
