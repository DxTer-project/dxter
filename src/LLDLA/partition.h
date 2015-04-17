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

#include "DLANode.h"
#include "node.h"
#include "layers.h"
#include "DLAOp.h"

#if DOLLDLA

class Partition : public DLANode
{
  
 private:
  Dir m_partType;

  unsigned int m_splitSize;

  const SizeList* m_startSizes;
  const SizeList* m_endSizes;

  DataTypeInfo* m_startInfo;
  DataTypeInfo* m_endInfo;

  Name m_startName;
  Name m_endName;

  Name m_startDim;
  Name m_endDim;

  void SetHorizontalNames();
  void SetVerticalNames();

  void BuildHorizontalDataTypeCache();
  void BuildVerticalDataTypeCache();

  void BuildHorizontalDataTypeInfo();
  void BuildVerticalDataTypeInfo();

  void BuildHorizontalSizes();
  void BuildVerticalSizes();
  void BuildStartAndEndSizes(const SizeList* toSplit);

 public:
  Layer m_layer;
  
  Partition(Layer layer, Dir partType, Size splitSize);

  virtual void PrintCode(IndStream &out);

  virtual ~Partition() {}
  inline void SetLayer(Layer layer) { m_layer = layer; }
  inline Layer GetLayer() const { return m_layer; }
  virtual bool IsReadOnly() const { return false; }
  virtual bool CanTrans() const { return false; }

  virtual NodeType GetType() const { return "Partition" + std::to_string((long long int) m_partType) + std::to_string((long long int) m_splitSize); }
  static ClassType GetClass() {return "partitionNode";}
  virtual ClassType GetNodeClass() const { return GetClass(); }
  virtual Name GetName(ConnNum num) const;

  virtual void AddVariables(VarSet &set) const;
  virtual void BuildDataTypeCache();
  virtual void ClearDataTypeCache();

  virtual ConnNum NumOutputs() const;
  virtual const DataTypeInfo& DataType(ConnNum num) const;

  // TODO LOOK FOR MULTIPLE PARTITIONS
  virtual bool Overwrites(const Node* input, ConnNum num) const { return false; }

  virtual Node* GetNewInst() { return BlankInst(); }
  static Node* BlankInst() { return new Partition(ABSLAYER, VERTICAL, 5); }

  virtual void Duplicate(const Node* orig, bool shallow, bool possMerging);

  virtual void Prop();

  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
};

class PartitionLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
  Size m_bs;

  PartitionLowerLayer(Layer fromLayer, Layer toLayer);
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

#endif //DOLLDLA
