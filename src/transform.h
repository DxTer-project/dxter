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

#include "base.h"
#include "node.h"
#include "universe.h"
#include "poss.h"
#include "pset.h"

//Main Transformation class
class Transformation
{
 public:
  virtual ~Transformation() {}
  virtual string GetType() const {return "Transformation";}
  virtual bool IsSingle() const {return false;}
  virtual bool IsRef() const {return false;}
  virtual bool IsVarRef() const {return false;}
};

//Holds one transformation
class SingleTrans : public Transformation
{
 public:
  SingleTrans() {}
  ~SingleTrans() {}
  virtual string GetType() const {return "SingleTrans";}
  virtual bool IsSingle() const {return true;}
  virtual bool CanApply(const Poss *poss, const Node *node) const = 0;
  virtual void Apply(Poss *poss, Node *node) const = 0;
  virtual bool WorthApplying(const Node *node) const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const {throw;}
};

//Holds one transformation
class VarTrans : public Transformation
{
 public:
  VarTrans() {}
  ~VarTrans() {}
  virtual string GetType() const {return "VarTrans";}
  virtual bool IsVarRef() const {return true;}
  virtual bool IsMultiRef() const {return false;}
  virtual int CanApply(const Poss *poss, const Node *node, void **cache) const = 0;
  virtual void Apply(Poss *poss, int num, Node *node, void **cache) const = 0;
  virtual bool WorthApplying(const Node *node) const {return true;}
  virtual Cost RHSCostEstimate(const Node *node) const {throw;}
  virtual void CleanCache(void **cache) const = 0;
};

//Holds multiple transformation and only allows the (estimated)
// best MAXNUMBEROFREFINEMENTS transformations
class MultiTrans : public VarTrans
{
 public:
  TransConstVec m_trans;
  bool m_isRef;
  ~MultiTrans() {}
  void AddTrans(SingleTrans *trans);
  virtual bool IsMultiRef() const {return true;}
  virtual string GetType() const {return "MultiTrans";}
  virtual TransConstVec* GetApplicableTrans(const Poss *poss, const Node *node) const;
  virtual bool IsRef() const {return m_isRef;}
  unsigned int NumTransformations() const {return m_trans.size();}
  virtual void CleanCache(void **cache) const;
  virtual int CanApply(const Poss *poss, const Node *node, void **cache) const;
  virtual void Apply(Poss *poss, int num, Node *node, void **cache) const;
  virtual const Transformation* GetTrans(void **cache, int num) const;
};

