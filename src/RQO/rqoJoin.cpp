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


#include "rqoJoin.h"


#if DORQO

Join::Join()
  : Sortable()
{

}

Join::Join(string sortBy, 
	   vector<string> in0Fields, 
	   vector<string> in1Fields)
  : Sortable(sortBy),
    m_in0Fields(in0Fields),
    m_in1Fields(in1Fields)
{
  static int num = 1;
  m_name = "join" + std::to_string(num);
  ++num;
}

NodeType Join::GetType() const
{
  string ret = m_sortBy;
  if (m_in0Fields.size() != m_in1Fields.size())
    throw;
  vector<string>::const_iterator iter0 = m_in0Fields.begin();
  vector<string>::const_iterator iter1 = m_in1Fields.begin();  
  for(; iter0 != m_in0Fields.end(); ++iter0, ++iter1) {
    ret += "," + *iter0 + *iter1;
  }
  return ret;
}

void Join::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  const Join *join = (Join*)orig;
  m_name = join->m_name;
  m_sortBy = join->m_sortBy;
  m_in0Fields = join->m_in0Fields;
  m_in1Fields = join->m_in1Fields;
  Node::Duplicate(orig, shallow, possMerging);
}

const DataTypeInfo& Join::DataType(ConnNum num) const
{
  return m_dataTypeInfo;
}

void Join::ClearDataTypeCache()
{
  
}

void Join::BuildDataTypeCache()
{
  m_dataTypeInfo.m_sortedBy = m_sortBy;
  m_dataTypeInfo.m_fields = InputDataType(0).m_fields;
  m_dataTypeInfo.m_fields.insert(InputDataType(1).m_fields.begin(),
				 InputDataType(1).m_fields.end());
}

void Join::Prop()
{
  //Node::Prop();
  if (m_inputs.size() != 2)
    throw;
  if (m_in0Fields.size() != m_in1Fields.size())
    throw;
  if (m_in0Fields.empty())
    throw;
  if (m_name.empty())
    throw;

  const DataTypeInfo &in0 = InputDataType(0);
  const DataTypeInfo &in1 = InputDataType(1);


  for (auto str : m_in0Fields) {
    if (in0.m_fields.find(str) == in0.m_fields.end()) {
      cout << "input A does not include field over which we are joining\n";
      throw;
    }
  }

  for (auto str : m_in1Fields) {
    if (in1.m_fields.find(str) == in1.m_fields.end()){
      cout << "input B does not include field over which we are joining\n";
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

void Join::PrintCode(IndStream &out)
{
  out.Indent();
  string in0 = GetInputNameStr(0);
  string in1 = GetInputNameStr(1);
  *out << m_name << " = BaseJoin( " << in0
       << ", " << in1;
  vector<string>::iterator iter0 = m_in0Fields.begin();
  vector<string>::iterator iter1 = m_in1Fields.begin();  
  for(; iter0 != m_in0Fields.end(); ++iter0, ++iter1) {
    *out << ", " << in0 << "." << *iter0 << " = "
	 << in1 << "." << *iter1;
  }
  *out << ", " << m_sortBy << " );\n";
}

Join* Join::CreateCopyOfJoin() const
{
  Join *newJoin = new Join(m_sortBy,
			   m_in0Fields,
			   m_in1Fields);
  return newJoin;
}

SwapNodes::SwapNodes(unsigned int inNum, ClassType type) 
  : m_inNum(inNum),
    m_type(type)
{
  if (inNum > 1)
    throw;
}

bool SwapNodes::CanApply(const Node *node) const
{
  throw;
  if (!node->IsJoin())
    throw;

  Join *join = (Join*)node;
  std::vector<string>::iterator it;
  Node *inNode = node->Input(m_inNum);
  //inputs to both joins have to match up

  if (inNode->IsJoin())
  {
    Join *inJoin = (Join*)inNode;
    if(m_inNum == 0)
    {
      for (auto str : join->m_in0Fields) 
      {
        for(it = inJoin->m_in0Fields.begin(); it != inJoin->m_in0Fields.end(); ++it)
        {
          if(*it == str)
          {
            break;
          }
        }
        if(it == inJoin->m_in0Fields.end())
        {
          return false;
        }
      }
    }
    else if(m_inNum == 1)
    {
      for (auto str : join->m_in1Fields) 
      {
        for(it = inJoin->m_in1Fields.begin(); it != inJoin->m_in1Fields.end(); ++it)
        {
          if(*it == str)
          {
            break;
          }
        }
        if(it == inJoin->m_in1Fields.end())
        {
          return false;
        }
      }
    }
    return true;
  }
  else
    return false;
}

void SwapNodes::Apply(Node *node) const
{
  throw;
  Join *inputJoin = (Join*)(node->Input(m_inNum));
  Join *newInputJoin = inputJoin->CreateCopyOfJoin();

  cout << "applying " << GetType() << " to poss with ";
  node->m_poss->PrintTransVecUp();

  node->m_poss->AddNode(newInputJoin);

  node->RedirectChildren(newInputJoin);
  if (m_inNum == 0)
    newInputJoin->AddInput(node, 0);
  else
    newInputJoin->AddInput(inputJoin->Input(0), inputJoin->InputConnNum(0));

  if (m_inNum == 1)
    newInputJoin->AddInput(node, 0);
  else
    newInputJoin->AddInput(inputJoin->Input(1), inputJoin->InputConnNum(1));
  
  
  node->ChangeInput2Way(inputJoin, 0,
			inputJoin->Input(m_inNum), inputJoin->InputConnNum(m_inNum));
  
  if (inputJoin->m_children.empty()) {
    inputJoin->m_poss->DeleteChildAndCleanUp(inputJoin);
  }
}

#endif
