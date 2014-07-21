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

#include "transform.h"
#include "DLAOp.h"


#if DOSQM || DOSM
#define TWOSIDEDTRXMCOMPONENTSLAYER ABSLAYER
#define TWOSIDEDTRXMLAYER S3LAYER
#elif DOELEM
#define TWOSIDEDTRXMCOMPONENTSLAYER DMLAYER
#define TWOSIDEDTRXMLAYER DMLAYER
#elif DOTENSORS||DOLLDLA
#define SKIPTWOSIDED 
#else
laksjdflkajsdf
#endif


#ifndef SKIPTWOSIDED
Loop* TwoSidedTrsmLowerVar1Alg(Node *in0, ConnNum num0,
			     Node *in1, ConnNum num1,
			    Layer layerBLAS, Layer layerTwoSidedTrxm);

Loop* TwoSidedTrsmLowerVar2Alg(Node *in0, ConnNum num0,
			     Node *in1, ConnNum num1,
			    Layer layerBLAS, Layer layerTwoSidedTrxm);

Loop* TwoSidedTrsmLowerVar4Alg(Node *in0, ConnNum num0,
			     Node *in1, ConnNum num1,
			     Layer layerBLAS, Layer layerTwoSidedTrxm);

Loop* TwoSidedTrmmLowerVar1Alg(Node *in0, ConnNum num0,
			     Node *in1, ConnNum num1,
			    Layer layerBLAS, Layer layerTwoSidedTrxm);

Loop* TwoSidedTrmmLowerVar2Alg(Node *in0, ConnNum num0,
			     Node *in1, ConnNum num1,
			    Layer layerBLAS, Layer layerTwoSidedTrxm);

Loop* TwoSidedTrmmLowerVar4Alg(Node *in0, ConnNum num0,
			    Node *in1, ConnNum num1,
			    Layer layerBLAS, Layer layerTwoSidedTrxm);


//Old name for two-sided Trsm/Trmm problem
//Left Hegst = two-sided Trmm
//Right Hegst = two-sided Trsm
class TwoSidedTrxm : public DLAOp<2,1>
{
 public:
  bool m_invert;
  Tri m_tri;
 TwoSidedTrxm(Layer layer, bool invert, Tri tri) 
   : m_invert(invert), m_tri(tri) {SetLayer(layer);}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "TwoSidedTrxm";}
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual NodeType GetType() const;
  static Node* BlankInst() { return  new TwoSidedTrxm(ABSLAYER, false, LOWER); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual Phase MaxPhase() const;
  virtual bool ShouldCullDP() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
};


#if DODPPHASE
class DistTwoSidedTrxmToLocalTwoSidedTrxm : public SingleTrans
{
 public:
  virtual string GetType() const {return "Distributed TwoSidedTrxm to Local TwoSidedTrxm";}
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};
#endif


class TwoSidedTrxmLoopExp : public SingleTrans
{
 public:
  int m_varNum;
  Layer m_fromLayer, m_toLayerBLAS, m_toLayerTwoSidedTrxm;
 TwoSidedTrxmLoopExp(int varNum, Layer fromLayer, Layer toLayerBLAS, Layer toLayerTwoSidedTrxm) 
   : m_varNum(varNum), m_fromLayer(fromLayer), 
    m_toLayerBLAS(toLayerBLAS), m_toLayerTwoSidedTrxm(toLayerTwoSidedTrxm) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};

class TwoSidedTrxmLowerLayer : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toLayer;
 TwoSidedTrxmLowerLayer(Layer fromLayer, Layer toLayer)
   : m_fromLayer(fromLayer), m_toLayer(toLayer) {}
  virtual string GetType() const { return "TwoSidedTrxm lower layer"; }
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
};

#endif // !SKIPTWOSIDED
