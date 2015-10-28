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



#pragma once

#include "layers.h"
#include "sortable.h"
#include "rqoBasis.h"
#include "transform.h"
#include <math.h>

#if DORQO

class Sort : public Sortable
{

 protected:
  DataTypeInfo m_dataTypeInfo;


 public:
  Sort();
  Sort(string sortBy);
  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new Sort; }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual const DataTypeInfo& DataType(ConnNum num) const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Cost GetCost() {return Input(0)->Outputs() * log(Input(0)->Outputs());}
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "sort";}
  virtual void ClearDataTypeCache();
  virtual void BuildDataTypeCache();
  virtual bool Overwrites(const Node *input, ConnNum num) const {return false;}
  virtual int Outputs() {return Input(0)->Outputs();}
};

//This transformation removes a sort node that is before a sort node.
class RemoveExtraSort : public SingleTrans
{
 public:
  virtual string GetType() const {return "Remove Extra Sort Nodes";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

//This transformation deletes a sortby in front of the node it is checking, since
//it would be redundant.
class RemoveRedundantSortBy : public SingleTrans
{
 public:
  virtual string GetType() const {return "Remove Redundant SortBy Fields";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

//This removes a node that is a sort before a sortable node. This is not
//superfluous since we allow sortable nodes to have null sortby fields.
class RemoveSortBeforeSortable : public SingleTrans
{
 public:
  virtual string GetType() const {return "Remove A Node Before A SortBy";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};




#endif
