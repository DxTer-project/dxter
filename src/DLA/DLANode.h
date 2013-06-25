/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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

//Ancestor class for all DLA nodes
class DLANode : public Node
{
 public:
  Cost m_cost;
  Layer m_layer;

  //Implement in sub-classes
  virtual DistType GetDistType(unsigned int num) const { return D_MC_MR; }
  virtual void SanityCheck();
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  
  DLANode();
  DLANode(Layer layer);
  inline void SetLayer(Layer layer) {m_layer = layer;}
  inline Layer GetLayer() const {return m_layer;}
  virtual bool IsReadOnly() const {return false;}
  virtual bool CanTrans() const {return false;}
  virtual const Sizes* GetM(unsigned int num) const = 0;
  virtual const Sizes* GetN(unsigned int num) const = 0;
  virtual const Sizes* LocalM(unsigned int num) const = 0;
  virtual const Sizes* LocalN(unsigned int num) const = 0;
  const Sizes* GetInputM(unsigned int num) const;
  const Sizes* GetInputN(unsigned int num) const;
  const Sizes* InputLocalM(unsigned int num) const;
  const Sizes* InputLocalN(unsigned int num) const;
  virtual void ClearBeforeProp();
  DistType InputDistType(unsigned int num) const;
  virtual string GetCostStr();
  virtual bool HasProped() const;
  virtual bool IsDLA() const {return true;}
  virtual Cost GetCost() {return m_cost;}
  DLANode* FindNonRedistParent(unsigned int num);
  DLANode* FindNonRedistParent(unsigned int num, unsigned int &parentNum);
  DLANode* FindSideEffectingUser(unsigned int num);
  bool IsRowVec(unsigned int num) const;
  bool IsColVec(unsigned int num) const;
  bool IsVec(unsigned int num) const;
  bool IsScalar(unsigned int num) const;
  virtual bool ShouldCullDP() const {return false;}
  virtual bool DoNotCullDP() const {return false;}
  virtual bool ShouldCullSR() const {return false;}
  virtual bool DoNotCullSR() const {return false;}
  virtual bool CanTransposeInputs() const {return false;} 
  virtual void UpdateInnerPackingMultiple(PackSize size);
  virtual bool IsBLISParallelizable() const {return false;}
  virtual void Parallelize(Comm comm) {throw;}
};

void DLACullDP(Poss *poss, bool &cullIfPossible, bool &doNotCull);
void DLACullRO(Poss *poss, bool &cullIfPossible, bool &doNotCull);
void DLACullLA(Poss *poss, bool &cullIfPossible, bool &doNotCull);
#if DOSQM || DOSM
void DLACullSR(Poss *poss, bool &cullIfPossible, bool &doNotCull);
#endif
