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
#include "rqoScan.h"
#include "rqoIndexedNode.h"
#include "rqoNIndexedNode.h"
#include "rqoOrderedIndexNode.h"

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
  m_query(query),
  m_fields(fields),
  m_sortBy(sortBy)
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
  m_relation = node->m_relation;
  m_fields = node->m_fields;
  m_sortBy = node->m_sortBy;
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
}

void InputNode::PrintCode(IndStream &out)
{
  out.Indent();
  *out << m_varName << " = scanFunc(" << m_fileName << "," << m_sortBy << "," << m_query << ");\n";
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


bool InputToScan::CanApply(const Node *node) const
{
  if(node->GetNodeClass() != "inputNode")
  {
    cout << "Throwing in InputToScan CanApply since : " << node->GetClass() << endl;
    throw;
  }

  return true;
}

void InputToScan::Apply(Node *node) const
{
  InputNode *orig = (InputNode*) node;
  Scan *newScan = new Scan(orig->m_varName, orig->m_sortBy, orig->m_fields, 
    orig->m_fileName, orig->m_query);
  newScan->SetRelation(orig->m_relation);

  orig->m_poss->AddNode(newScan);
  orig->RedirectChildren(newScan);
  if(orig->m_children.empty())
  {
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
}

bool InputToIndex::CanApply(const Node *node) const
{
  if(node->GetNodeClass() != "inputNode")
  {
    cout << "Throwing in InputToIndex CanApply since : " << node->GetClass() << endl;
    throw;
  }

  InputNode *newNode = (InputNode*)node;
  if(newNode->m_relation->indeces.size() != 1)
  {
    throw;
  }

  return true;
}

void InputToIndex::Apply(Node *node) const
{
  InputNode *orig = (InputNode*) node;
  set<int>::iterator iter = orig->m_relation->indeces.begin();
  IndexedNode *newIndex = new IndexedNode(orig->m_varName, orig->m_sortBy, orig->m_fields, 
    orig->m_fileName, orig->m_query, (*iter));
  newIndex->SetRelation(orig->m_relation);

  orig->m_poss->AddNode(newIndex);
  orig->RedirectChildren(newIndex);
  if(orig->m_children.empty())
  {
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
}

bool InputToNIndex::CanApply(const Node *node) const
{
  if(node->GetNodeClass() != "inputNode")
  {
    cout << "Throwing in InputToNIndex CanApply since : " << node->GetClass() << endl;
    throw;
  }

  InputNode *newNode = (InputNode*)node;
  if(newNode->m_relation->indeces.size() <= 1)
  {
    throw;
  }

  return true;
}

void InputToNIndex::Apply(Node *node) const
{
  InputNode *orig = (InputNode*) node;
  NIndexedNode *newIndex = new NIndexedNode(orig->m_varName, orig->m_sortBy, orig->m_fields, 
    orig->m_fileName, orig->m_query, orig->m_relation->indeces);
  newIndex->SetRelation(orig->m_relation);

  orig->m_poss->AddNode(newIndex);
  orig->RedirectChildren(newIndex);
  if(orig->m_children.empty())
  {
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
}

bool InputToOrdered::CanApply(const Node *node) const
{
  if(node->GetNodeClass() != "inputNode")
  {
    cout << "Throwing in InputToOrdered CanApply since : " << node->GetClass() << endl;
    throw;
  }

  InputNode *newNode = (InputNode*)node;
  if(newNode->m_relation->indeces.size() != 1)
  {
    throw;
  }

  return true;
}

void InputToOrdered::Apply(Node *node) const
{
  InputNode *orig = (InputNode*) node;
  set<int>::iterator iter = orig->m_relation->indeces.begin();
  OrderedIndexedNode *newIndex = new OrderedIndexedNode(orig->m_varName, orig->m_sortBy, orig->m_fields, 
    orig->m_fileName, orig->m_query, (*iter));
  newIndex->SetRelation(orig->m_relation);

  orig->m_poss->AddNode(newIndex);
  orig->RedirectChildren(newIndex);
  if(orig->m_children.empty())
  {
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
}


#endif //DORQO

