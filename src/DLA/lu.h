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
#if DOBLIS||DOELEM

#include "base.h"
#include "DLAOp.h"
#include "DLANode.h"

RealLoop* LUVar5Loop(Node *Ain, ConnNum Anum,
		 Node *Pin, ConnNum Pnum,
		 Layer BLASLayer, Layer LAPACKLayer);

class LU : public DLAOp<2,2>
{
 public:
  LU(Layer layer) {SetLayer(layer);}
  static Node* BlankInst() { return new LU(ABSLAYER); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "LU";}
  virtual NodeType GetType() const {return "LU";}
  virtual Phase MaxPhase() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
};


class PanelLU : public DLAOp<7,7>
{
 public:
  PanelLU(Layer layer) {SetLayer(layer);}
  static Node* BlankInst() { return new PanelLU(ABSLAYER); }
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "PanelLU";}
  virtual NodeType GetType() const {return "PanelLU";}
  virtual Phase MaxPhase() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
};

class LULoopExp : public SingleTrans
{
 public:
  Layer m_fromLayer, m_toBLASLayer, m_toLAPACKLayer;
  unsigned int m_var;
 LULoopExp(Layer fromLayer, Layer toBLASLayer, Layer toLAPACKLayer, unsigned int var) 
   : m_fromLayer(fromLayer), m_toBLASLayer(toBLASLayer), 
    m_toLAPACKLayer(toLAPACKLayer), m_var(var) {}
  virtual string GetType() const;
  virtual bool CanApply(const Node *node) const;
  virtual void Apply(Node *node) const;
  virtual bool IsRef() const {return true;}
};
#endif
