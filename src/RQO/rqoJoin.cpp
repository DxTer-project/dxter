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
#include "hJoin.h"
#include "nJoin.h"
#include "mJoin.h"


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
  string ret = GetClass() + " " + m_sortBy;
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
      cout << "input A does not include field over which we are joining:" << str << endl;
      cout << "m_fields.size() = " << in0.m_fields.size() << endl;
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
  *out << m_name << " = nestedJoin(" << m_sortBy << ","
    << in0 << "," << in1;
  vector<string>::iterator iter0 = m_in0Fields.begin();
  vector<string>::iterator iter1 = m_in1Fields.begin();  
  for(; iter0 != m_in0Fields.end(); ++iter0, ++iter1) {
    *out << "," << in0 << "." << *iter0 << ","
   << in1 << "." << *iter1;
  }
  *out << ");\n";
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
  if (!node->IsJoin())
    throw;

  Join *join = (Join*)node;
  std::vector<string>::iterator it;
  Node *inNode = node->Input(m_inNum);
  //we start with swappable as true, and will attempt to prove it wrong
  bool swappable = true;


  if (inNode->IsJoin())
  {
    Join *inJoin = (Join*)inNode;
    //collect the inputs from the join we want to swap with
    const DataTypeInfo &in0 = inJoin->InputDataType(0);
    const DataTypeInfo &in1 = inJoin->InputDataType(1);
    //swapping on the root nodes first input. the contents of m_inNum == 1 is 
    //identical to this one
    if(m_inNum == 0)
    {
    	//The goal here is to find one of the inputs that matches everything
    	//we need entirely. If the first one fails, we set swappable to false
    	//and try the other one.

      //we have to make sure that every field that was part of
      // this join is coming from the same place, so first check in0 of input
      for (auto str : join->m_in0Fields) 
      {
        if (in0.m_fields.find(str) == in0.m_fields.end()) {
          swappable = false;
          break;
        }
      }
      //this if makes sure that if the sortby doesnt exist where we're swapping
      //we dont try it
      if (!join->m_sortBy.empty() 
	  && (in0.m_fields.find(join->m_sortBy) == in0.m_fields.end())
	  && (join->InputDataType(1).m_fields.find(join->m_sortBy) == join->InputDataType(1).m_fields.end()))
	{
	  swappable = false;
	}

      //if swappable is false, we try the second input
      //otherwise, tell the node what to swap with
      if(!swappable)
      {
        swappable = true;
        for(auto str : join->m_in0Fields)
        {
          if (in1.m_fields.find(str) == in1.m_fields.end()) {
            swappable = false;
            break;
          }
        }
	if (!join->m_sortBy.empty() 
	    && (in1.m_fields.find(join->m_sortBy) == in1.m_fields.end())
	    && (join->InputDataType(1).m_fields.find(join->m_sortBy) == join->InputDataType(1).m_fields.end()) )
	  {
	    swappable = false;
	  }
      }

      return swappable;
    }
    else if(m_inNum == 1)
    {
      for (auto str : join->m_in1Fields) 
      {
        if (in0.m_fields.find(str) == in0.m_fields.end()) {
          swappable = false;
          break;
        }
      }
      
      if (!join->m_sortBy.empty() 
	  && (in0.m_fields.find(join->m_sortBy) == in0.m_fields.end())
	  && (join->InputDataType(0).m_fields.find(join->m_sortBy) == join->InputDataType(0).m_fields.end()))
	{
	  swappable = false;
	}

      //if swappable is false, we try the second input
      //otherwise, tell the node what to swap with
      if(!swappable)
      {
        swappable = true;
        for(auto str : join->m_in1Fields)
        {
          if (in1.m_fields.find(str) == in1.m_fields.end()) {
            swappable = false;
            break;
          }
        }

	if (!join->m_sortBy.empty() 
	    && (in1.m_fields.find(join->m_sortBy) == in1.m_fields.end())
	    && (join->InputDataType(0).m_fields.find(join->m_sortBy) == join->InputDataType(0).m_fields.end()))
	  {
	    swappable = false;
	  }

      }

      return swappable;
    }
  }
  return false;
}

void SwapNodes::Apply(Node *node) const
{
  Join *inputJoin = (Join*)(node->Input(m_inNum));
  Join *newInputJoin = inputJoin->CreateCopyOfJoin();

  //cout << "applying " << GetType() << " to poss with ";
  //node->m_poss->PrintTransVecUp();

  node->m_poss->AddNode(newInputJoin);

  node->RedirectChildren(newInputJoin);

  //recalculate everything from canapply
  Join *join = (Join*)node;
  //input to keep is the input that we want to keep going into root node when we're finished
  int inputToKeep;
  bool swappable = true;
  const DataTypeInfo &in0 = inputJoin->InputDataType(0);
  const DataTypeInfo &in1 = inputJoin->InputDataType(1);
  //swapping on input 0
  //everything in these if statements are doing everythign we did from can apply, 
  //but this time with the goal of saving which input we want to keep
  //Once again, everything in m_inNum == 1 is identical to this just for the other input
  if(m_inNum == 0)
    {
      for (auto str : join->m_in0Fields) 
	{
	  if (in0.m_fields.find(str) == in0.m_fields.end()) {
	    swappable = false;
	    break;
	  }
	}
      if(swappable)
	{
	  inputToKeep = 0;
	}
      else
	{
	  swappable = true;
	  for(auto str : join->m_in0Fields)
	    {
	      if (in1.m_fields.find(str) == in1.m_fields.end()) {
		swappable = false;
		break;
	      }
	    }
	  if(swappable)
	    {
	      inputToKeep = 1;
	    }
	  else
	    throw;
	}
	//This code is where we put in new inputs to each of the nodes we are swapping
      if(inputToKeep == 0)
	{
	  newInputJoin->AddInput(node, 0);
	  newInputJoin->AddInput(inputJoin->Input(1), inputJoin->InputConnNum(1));

	  node->ChangeInput2Way(inputJoin, 0,
				inputJoin->Input(inputToKeep), inputJoin->InputConnNum(inputToKeep));
	}
      else
	{
	  newInputJoin->AddInput(node, 0);
	  newInputJoin->AddInput(inputJoin->Input(0), inputJoin->InputConnNum(0));


	  node->ChangeInput2Way(inputJoin, 0,
				inputJoin->Input(inputToKeep), inputJoin->InputConnNum(inputToKeep));
	}
    }
  else if(m_inNum == 1)
    {
      for (auto str : join->m_in1Fields) 
	{
	  if (in0.m_fields.find(str) == in0.m_fields.end()) {
	    swappable = false;
	    break;
	  }
	}
      //if swappable is false, we try the second input
      //otherwise, tell the node what to swap with
      if(swappable)
	{
	  inputToKeep = 0;
	}
      else
	{
	  swappable = true;
	  for(auto str : join->m_in1Fields)
	    {
	      if (in1.m_fields.find(str) == in1.m_fields.end()) {
		swappable = false;
		break;
	      }
	    }
	  if(swappable)
	    {
	      inputToKeep = 1;
	    }
	  else
	    throw;
	}


      if(inputToKeep == 0)
	{
	  newInputJoin->AddInput(inputJoin->Input(1), inputJoin->InputConnNum(1));
	  newInputJoin->AddInput(node, 0);

	  node->ChangeInput2Way(inputJoin, 0,
				inputJoin->Input(inputToKeep), inputJoin->InputConnNum(inputToKeep));
	}
      else
	{
	  newInputJoin->AddInput(inputJoin->Input(0), inputJoin->InputConnNum(0));
	  newInputJoin->AddInput(node, 0);

	  node->ChangeInput2Way(inputJoin, 0,
				inputJoin->Input(inputToKeep), inputJoin->InputConnNum(inputToKeep));
	}
    }
  

  
  if (inputJoin->m_children.empty()) {
    inputJoin->m_poss->DeleteChildAndCleanUp(inputJoin);
  }
}


bool JoinToHash::CanApply(const Node *node) const
{
  if(node->GetNodeClass() != "join")
  {
    cout << "Throwing in JoinToHash CanApply since : " << node->GetClass() << endl;
    throw;
  }

  return true;
}

void JoinToHash::Apply(Node *node) const
{
  Join *orig = (Join*) node;
  HJoin *newJoin = new HJoin(orig->m_sortBy, orig->m_in0Fields, orig->m_in1Fields);

  orig->m_poss->AddNode(newJoin);
  orig->RedirectChildren(newJoin);
  newJoin->AddInput(orig->Input(0), orig->InputConnNum(0));
  newJoin->AddInput(orig->Input(1), orig->InputConnNum(1));
  if(orig->m_children.empty())
  {
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
}

bool JoinToNested::CanApply(const Node *node) const
{
  if(node->GetNodeClass() != "join")
    throw;

  return true;
}

void JoinToNested::Apply(Node *node) const
{
  Join *orig = (Join*) node;
  NJoin *newJoin = new NJoin(orig->m_sortBy, orig->m_in0Fields, orig->m_in1Fields);

  orig->m_poss->AddNode(newJoin);
  orig->RedirectChildren(newJoin);
  newJoin->AddInput(orig->Input(0), orig->InputConnNum(0));
  newJoin->AddInput(orig->Input(1), orig->InputConnNum(1));
  if(orig->m_children.empty())
  {
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
}

bool JoinToMerge::CanApply(const Node *node) const
{
  if(node->GetNodeClass() != "join")
    throw;

  return true;
}

void JoinToMerge::Apply(Node *node) const
{
  Join *orig = (Join*) node;
  MJoin *newJoin = new MJoin(orig->m_sortBy, orig->m_in0Fields, orig->m_in1Fields);

  orig->m_poss->AddNode(newJoin);
  orig->RedirectChildren(newJoin);
  newJoin->AddInput(orig->Input(0), orig->InputConnNum(0));
  newJoin->AddInput(orig->Input(1), orig->InputConnNum(1));
  if(orig->m_children.empty())
  {
    orig->m_poss->DeleteChildAndCleanUp(orig);
  }
}

#endif
