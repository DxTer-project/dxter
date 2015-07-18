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


#include "crossProduct.h"


#if DORQO

CrossProduct::CrossProduct()
    : RQONode()
    {
        static int num = 1;
        m_name = "crossproduct" + std::to_string(num);
        ++num;
    }

/*CrossProduct::CrossProduct(vector<string> in0Fields, 
       vector<string> in1Fields)
  : m_in0Fields(in0Fields),
    m_in1Fields(in1Fields)
{
  static int num = 1;
  m_name = "crossproduct" + std::to_string(num);
  ++num;
}*/

NodeType CrossProduct::GetType() const
{
  string ret = m_name;
  /*if (m_in0Fields.size() != m_in1Fields.size())
    throw;
  vector<string>::const_iterator iter0 = m_in0Fields.begin();
  vector<string>::const_iterator iter1 = m_in1Fields.begin();  
  for(; iter0 != m_in0Fields.end(); ++iter0, ++iter1) {
    ret += "," + *iter0 + *iter1;
  }*/
  return ret;
}

void CrossProduct::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const CrossProduct *cProd = (CrossProduct*)orig;
  m_name = cProd->m_name;
  //m_in0Fields = cProd->m_in0Fields;
  //m_in1Fields = cProd->m_in1Fields;
  Node::Duplicate(orig, shallow, possMerging);
}

const DataTypeInfo& CrossProduct::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void CrossProduct::ClearDataTypeCache()
{
  
}

void CrossProduct::BuildDataTypeCache()
{
  m_dataTypeInfo.m_sortedBy = "";
  m_dataTypeInfo.m_fields = InputDataType(0).m_fields;
  m_dataTypeInfo.m_fields.insert(InputDataType(1).m_fields.begin(),
                 InputDataType(1).m_fields.end());
}

void CrossProduct::Prop()
{
  //Node::Prop();
  if (m_inputs.size() != 2)
    throw;
  if (m_name.empty())
    throw;



  
}

void CrossProduct::PrintCode(IndStream &out)
{
  out.Indent();
  string in0 = GetInputNameStr(0);
  string in1 = GetInputNameStr(1);
  *out << m_name << " = CrossProduct( " << in0
       << ", " << in1;
  /*vector<string>::iterator iter0 = m_in0Fields.begin();
  vector<string>::iterator iter1 = m_in1Fields.begin();  
  for(; iter0 != m_in0Fields.end(); ++iter0, ++iter1) {
    *out << ", " << in0 << "." << *iter0 << " = "
     << in1 << "." << *iter1;
  }*/
  *out << ");\n";
}

#endif