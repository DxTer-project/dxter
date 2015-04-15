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

DataTypeInfo Not::m_info;

string Not::GetUniqueName()
{
  static int num = 0;
  return (string)"not" + std::to_string(num++);
}

Not::Not()
:m_name(GetUniqueName())
{

}

void Not::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  Node::Duplicate(orig, shallow, possMerging);
  const Not *notOrig = (Not*)orig;
  m_name = notOrig->m_name;
}

Name Not::GetName(ConnNum num) const
{
  return m_name;
}

void Not::Prop()
{
  if (m_inputs.size() != 1)
    throw;
}

void Not::PrintCode(IndStream &out)
{
  out.Indent();
  *out << m_name.str() << " = Not( " << GetInputNameStr(0) << " );\n";
}



bool NotTrue::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Not::GetClass())
    throw;
  const Not *notNode = (Not*)node;
  if (notNode->Input(0)->GetNodeClass() == True::GetClass())
    return true;
  else
    return false;
}

void NotTrue::Apply(Node *node) const
{
  False *falseVal = new False;
  node->RedirectChildren(falseVal, 0);
  node->m_poss->AddNode(falseVal);
  node->m_poss->DeleteChildAndCleanUp(node);
}


bool NotFalse::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != Not::GetClass())
    throw;
  const Not *notNode = (Not*)node;
  if (notNode->Input(0)->GetNodeClass() == False::GetClass())
    return true;
  else
    return false;
}

void NotFalse::Apply(Node *node) const
{
  True *trueVal = new True;
  node->RedirectChildren(trueVal, 0);
  node->m_poss->AddNode(trueVal);
  node->m_poss->DeleteChildAndCleanUp(node);
}

#endif // DOBOOL

