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

#include "layers.h"

#if DOTENSORS

#include "transform.h"
#include "DLAOp.h"
#include "tensorRedist.h"
#include "lowerLayer.h"

class ZAxpBy : public DLAOp<3,1>
{
 public:
  Coef m_alpha, m_beta;
  DimVec m_permutation;
  ZAxpBy(Layer layer, Coef alpha, Coef beta);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual NodeType GetType() const;
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "zaxpby";}
  static Node* BlankInst() { return new ZAxpBy(ABSLAYER, COEFZERO, COEFZERO); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Phase MaxPhase() const;
  //  virtual bool ShouldCullDP() const;
  virtual bool DoNotCullDP() const;
};


class ZAxpByLowerLayer : public LowerLayer
{
 public:
  ZAxpByLowerLayer(Layer fromLayer, Layer toLayer)
   : LowerLayer(fromLayer, toLayer) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};


#endif //DOTENSORS

