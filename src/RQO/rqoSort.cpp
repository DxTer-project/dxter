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

Sort::Sort()
  : Sortable()
{

}

Sort::Sort(string sortBy)
  : Sortable(sortBy)
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


  if(m_sortBy.empty()) 
  {
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
  if (!node->IsSortable())
    throw;

  Sortable *sorted = (Sortable*)node;
  const Node *inNode = sorted->Input(0);
  if (!inNode->IsSortable())
    return false;
  const Sortable *input = (Sortable*)inNode;
  if (!sorted->m_sortBy.empty() && !input->m_sortBy.empty()) {
    return true;
  }
  else
    return false;
}

void RemoveRedundantSortBy::Apply(Node *node) const
{
  Sortable *input = (Sortable*)node->Input(0);
  input->m_sortBy = "";
}

bool RemoveSortBeforeSortable::CanApply(const Node *node) const
{
  if (node->m_inputs.size() != 1)
    throw;
  if (!node->IsSortable())
    throw;

  Sortable *sorted = (Sortable*)node;
  const Node *inNode = sorted->Input(0);
  if (!sorted->m_sortBy.empty() 
      && inNode->GetNodeClass() == Sort::GetClass()) 
    {
      return true;
    }
  else
    return false;
}

void RemoveSortBeforeSortable::Apply(Node *node) const
{
  Node *input = node->Input(0);
  node->ChangeInput2Way(input, node->InputConnNum(0),
            input->Input(0), input->InputConnNum(0));
  if (input->m_children.empty()) {
    input->m_poss->DeleteChildAndCleanUp(input);
  }
}


template <class T>
SwapNodes<T>::SwapNodes(unsigned int inNum) 
  : m_inNum(inNum) 
{
  if (inNum > 1)
    throw;
}

template <class T>
bool SwapNodes<T>::CanApply(const Node *node) const
{
  if (!node->IsJoin())
    throw;

  throw; // Add additional checks

  if (node->Input(m_inNum)->IsJoin())
    return true;
  else
    return false;
}

template <class T>
void SwapNodes<T>::Apply(Node *node) const
{
  Join *inputJoin = (Join*)(node->Input(m_inNum));
  Join *newInputJoin = inputJoin->CreateCopyOfJoin();
  node->RedirectChildren(newInputJoin);
  if (m_inNum == 0)
    newInputJoin->AddInput(node, 0);
  else
    newInputJoin->AddInput(inputJoin->Input(0), inputJoin->InputConnNum(0));

  if (m_inNum == 1)
    newInputJoin->AddInput(node, 0);
  else
    newInputJoin->AddInput(inputJoin->Input(1), inputJoin->InputConnNum(1));


  throw; //check numbers  
  node->ChangInput2Way(inputJoin, 0,
			 inputJoin->Input(m_inNum), inputJoin->InputConnNum(m_inNum));

  if (inputJoin->m_children.empty()) {
    inputJoin->m_poss->DeleteChildAndCleanUp(inputJoin);
  }
}


#endif

