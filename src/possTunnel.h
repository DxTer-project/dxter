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

#include "base.h"
#include "DLANode.h"

class PSet;

class PossTunnel : public DLANode
{
 public:
  PossTunType m_tunType;
  PSet *m_pset;
  PossTunnel();
  PossTunnel(PossTunType type);
  void SetPSet(PSet *set);
  virtual DistType GetDistType(unsigned int num) const;
  static Node* BlankInst() { return new PossTunnel;}
  virtual Node* GetNewInst() {return BlankInst(); }
  virtual PossTunnel* GetSetTunnel();
  virtual void Prop();
  virtual void PrintCode(IndStream &out) {};
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual void SanityCheck();
  virtual bool IsPossTunnel() const {return true;}
  virtual bool IsPossTunnel(PossTunType type) const;
  virtual unsigned int NumOutputs() const;
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "PossTunnel";}
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  virtual Name GetName(unsigned int num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, unsigned int num) const;
  virtual bool KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const;
};

string TunTypeToStr(PossTunType type);
