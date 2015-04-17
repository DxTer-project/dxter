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

#include "bool.h"

#if DOBOOL

DataTypeInfo And::m_info;

string And::GetUniqueName()
{
  static int num = 0;
  return (string)"and" + std::to_string(num++);
}

And::And()
:m_name(GetUniqueName())
{

}

void And::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  Node::Duplicate(orig, shallow, possMerging);
  const And *andOrig = (And*)orig;
  m_name = andOrig->m_name;
}

Name And::GetName(ConnNum num) const
{
  return m_name;
}

void And::Prop()
{
  if (m_inputs.size() != 2)
    throw;
}

void And::PrintCode(IndStream &out)
{
  out.Indent();
  *out << m_name.str() << " = And( " << GetInputNameStr(0)
       << " , " << GetInputNameStr(1) << " );\n";
}

bool AndToOr::CanApply(const Node *node) const
{
  return node->GetNodeClass() == And::GetClass();
}

void AndToOr::Apply(Node *node) const
{
  Not *not1 = new Not;
  not1->AddInput(node->Input(0), node->InputConnNum(0));

  Not *not2 = new Not;
  not2->AddInput(node->Input(1), node->InputConnNum(1));

  Or *newOr = new Or;
  newOr->AddInputs0(2, not1, not2);
  
  Not *not3 = new Not;
  not3->AddInput(newOr, 0);

  node->RedirectChildren(not3,0);

  Poss *poss = node->m_poss;
  poss->AddNodes(4, not1, not2, newOr, not3);

  poss->DeleteChildAndCleanUp(node); 
}



#endif // DOBOOL

