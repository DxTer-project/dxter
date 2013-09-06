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

#define LIMITEDSANITYCHECK

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
  virtual TransConstVec GetApplicableTrans(const Poss *poss, const Node *node) const = 0;
  virtual bool IsSingle() const {return false;}
  virtual bool IsRef() const {return false;}
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
  virtual TransConstVec GetApplicableTrans(const Poss *poss, const Node *node) const;
  virtual Cost RHSCostEstimate(const Node *node) const {throw;}
};

//Holds multiple transformation and only allows the (estimated)
// best MAXNUMBEROFREFINEMENTS transformations
class MultiTrans : public Transformation
{
 public:
  TransConstVec m_trans;
  bool m_isRef;
  ~MultiTrans();
  void AddTrans(SingleTrans *trans);
  virtual string GetType() const {return "MultiTrans";}
  virtual TransConstVec GetApplicableTrans(const Poss *poss, const Node *node) const;
  virtual bool IsRef() const {return m_isRef;}
  unsigned int NumTransformations() const {return m_trans.size();}
};

