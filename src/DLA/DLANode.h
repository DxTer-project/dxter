/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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
  virtual const Sizes* GetM(unsigned int num) const = 0;
  virtual const Sizes* GetN(unsigned int num) const = 0;
#if DODM
  virtual const Sizes* LocalM(unsigned int num) const = 0;
  virtual const Sizes* LocalN(unsigned int num) const = 0;
#endif
  const Sizes* GetInputM(unsigned int num) const;
  const Sizes* GetInputN(unsigned int num) const;
#if DODM
  const Sizes* InputLocalM(unsigned int num) const;
  const Sizes* InputLocalN(unsigned int num) const;
#endif
#elif DOTENSORS
  virtual const Dim NumDims(unsigned int num) const = 0;
  virtual const Sizes* Len(unsigned int num, Dim dim) const = 0;
  virtual const Sizes* LocalLen(unsigned int num, Dim dim) const = 0;  
  virtual const Dim InputNumDims(unsigned int num) const;
  const Sizes* InputLen(unsigned int num, Dim dim) const;
  const Sizes* InputLocalLen(unsigned int num, Dim dim) const;
#endif
  virtual void ClearBeforeProp();
  virtual bool IsDLA() const {return true;}
  virtual Cost GetCost() {return m_cost;}
#if DOELEM
  DLANode* FindNonRedistParent(unsigned int num);
  DLANode* FindNonRedistParent(unsigned int num, unsigned int &parentNum);
  virtual bool CanTransposeInputs() const {return false;} 
  virtual bool ShouldCullDP() const {return false;}
  virtual bool DoNotCullDP() const {return false;}
#endif
  DLANode* FindSideEffectingUser(unsigned int num);
#if TWOD
  bool IsScalar(unsigned int num) const;
#endif
#if DOBLIS
  virtual void UpdateInnerPackingMultiple(PackSize size);
  virtual bool IsBLISParallelizable() const {return false;}
  virtual void Parallelize(Comm comm) {throw;}
#endif
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


#if DODM

class DataTypeInfo
{
 public:
  DistType m_dist;
  DataTypeInfo();
  DataTypeInfo(DistType dist);
  DataTypeInfo(const DataTypeInfo &rhs);
  DataTypeInfo& operator=(const DataTypeInfo &rhs);
};

#elif !DOLLDLA
class DataTypeInfo
{
 public:
  DataTypeInfo& operator=(const DataTypeInfo &rhs) {return *this;}
};

#endif
