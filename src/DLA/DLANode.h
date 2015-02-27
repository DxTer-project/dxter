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
#include "node.h"
#include "comm.h"
#include "layers.h"

//Ancestor class for all DLA nodes
class DLANode : public Node
{
 public:
  Cost m_cost;
  Layer m_layer;

  //Implement in sub-classes
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  
  DLANode();
  DLANode(Layer layer);
  virtual ~DLANode() {}
  inline void SetLayer(Layer layer) {m_layer = layer;}
  inline Layer GetLayer() const {return m_layer;}
  virtual bool IsReadOnly() const {return false;}
  virtual bool CanTrans() const {return false;}

  virtual void Prop();

  //BAMTODO: eventually use DataTypeInfo to pass around data sizes

#if TWOD
  virtual const Sizes* GetM(ConnNum num) const = 0;
  virtual const Sizes* GetN(ConnNum num) const = 0;
  Size MaxNumberOfElements(ConnNum num) const;
#if DODM
  virtual const Sizes* LocalM(ConnNum num) const = 0;
  virtual const Sizes* LocalN(ConnNum num) const = 0;
#endif
  const Sizes* GetInputM(ConnNum num) const;
  const Sizes* GetInputN(ConnNum num) const;
#if DODM
  const Sizes* InputLocalM(ConnNum num) const;
  const Sizes* InputLocalN(ConnNum num) const;
#endif
  bool IsInputRowVector(ConnNum num) const;
  bool IsInputColVector(ConnNum num) const;
  bool IsInputScalar(ConnNum num) const;

#elif DOTENSORS
  virtual const Dim NumDims(ConnNum num) const = 0;
  virtual const Sizes* Len(ConnNum num, Dim dim) const = 0;
  virtual const Sizes* LocalLen(ConnNum num, Dim dim) const = 0;
  virtual const Dim InputNumDims(ConnNum num) const;
  const Sizes* InputLen(ConnNum num, Dim dim) const;
  const Sizes* InputLocalLen(ConnNum num, Dim dim) const;
  Size TotalNumberOfLocalElements(ConnNum num) const;
  Size TotalNumberOfInputLocalElements(ConnNum num) const;
  Size TotalNumberOfElements(ConnNum num) const;
  Size MaxNumberOfLocalElements(ConnNum num) const;
#endif
  virtual void ClearBeforeProp();
  virtual bool IsDLA() const {return true;}
  virtual Cost GetCost() {return m_cost;}
#if DOELEM
  DLANode* FindNonRedistParent(ConnNum num);
  DLANode* FindNonRedistParent(ConnNum num, ConnNum &parentNum);
  virtual bool CanTransposeInputs() const {return false;} 
#endif
#if DODM
  virtual bool DoNotCullDP() const {return false;}
#endif
#if DOELEM
  virtual bool ShouldCullDP() const {return false;}
#elif DOTENSORS
  virtual bool ShouldCullDP() const {return GetLayer() == DM1LAYER || GetLayer() == DM2LAYER;}
#endif
  DLANode* FindSideEffectingUser(ConnNum num);
#if TWOD
  void CheckInputNum(ConnNum num) const;
  bool IsScalar(ConnNum num) const;
#endif
#if DOBLIS
  virtual void UpdateInnerPackingMultiple(PackSize size);
  virtual bool IsBLISParallelizable() const {return false;}
  virtual void Parallelize(Comm comm) {throw;}
#endif
#if DOTENSORS
  virtual void AlignInfo(string &align,
			 DimVec &alignModes,
			 DimVec &alignModesSrc) {}
#endif

#if DOLLDLA

  int GetInputNumCols(ConnNum num) const;
  int GetInputNumRows(ConnNum num) const;

  int GetInputRowStride(ConnNum num) const;
  int GetInputColStride(ConnNum num) const;

  bool InputIsContiguous(ConnNum num) const;
  bool InputsAreSameSize(ConnNum left, ConnNum right) const;

  bool InputNIsMultipleOfVecRegWidth(ConnNum num) const;
  bool InputMIsMultipleOfVecRegWidth(ConnNum num) const;

#endif // DOLLDLA
};

#if TWOD

#if DOELEM
void DLACullDP(Poss *poss, bool &cullIfPossible, bool &doNotCull);
void DLACullRO(Poss *poss, bool &cullIfPossible, bool &doNotCull);
#endif
void DLACullLA(Poss *poss, bool &cullIfPossible, bool &doNotCull);
#if DOSQM || DOSM
void DLACullSR(Poss *poss, bool &cullIfPossible, bool &doNotCull);
#endif
#if DOLLDLA
void LLDLACull(Poss *poss, bool &cullIfPossible, bool &doNotCull);
#endif
#elif DOTENSORS
void TenCullDP(Poss *poss, bool &cullIfPossible, bool &doNotCull);
void TenCullRO(Poss *poss, bool &cullIfPossible, bool &doNotCull);

#endif


#if DOELEM

class DataTypeInfo
{
 public:
  DistType m_dist;
  DataTypeInfo();
  DataTypeInfo(DistType dist);
  DataTypeInfo(const DataTypeInfo &rhs);
  DataTypeInfo& operator=(const DataTypeInfo &rhs);
  const DistType& GetDist() const {return m_dist;}
  bool operator==(const DataTypeInfo &rhs) const  {return m_dist == rhs.m_dist;}
  bool operator!=(const DataTypeInfo &rhs) const  {return !(*this == rhs);}
};

#elif DOTENSORS
class DataTypeInfo
{
  DistType m_dist;
  Permutation m_perm;
 public:
  DataTypeInfo();
  DataTypeInfo(DistType dist);
  DataTypeInfo(const DataTypeInfo &rhs);
  DataTypeInfo(DistType dist, const Permutation &perm);
  DataTypeInfo& operator=(const DataTypeInfo &rhs);
  bool operator==(const DataTypeInfo &rhs) const;
  bool operator!=(const DataTypeInfo &rhs) const;
  DistType GetEffectiveDist() const;
  inline const DistType& GetDist() const {return m_dist;}
  inline const Permutation& GetPerm() const {return m_perm;}
  inline const bool HasPerm() const {return m_perm.HasPerm();}
  void SetToDefault(Dim numDims);
  void SetDistAndClearPerm(const DistType &dist);
  void SetPerm(const Permutation &perm);
  string Str() const;
};

#elif !DOLLDLA
class DataTypeInfo
{
 public:
  DataTypeInfo& operator=(const DataTypeInfo &rhs) {return *this;}
};

#endif
