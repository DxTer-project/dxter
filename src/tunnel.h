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

#include "base.h"
#include "DLANode.h"

class BasePSet;

class Tunnel : public DLANode
{
 public:
  TunType m_tunType;
  BasePSet *m_pset;
  Tunnel();
  Tunnel(TunType type);
  void SetPSet(BasePSet *set);
  static Node* BlankInst() { return new Tunnel;}
  virtual Node* GetNewInst() {return BlankInst(); }
  virtual Tunnel* GetSetTunnel();
  virtual void Prop();
  virtual void PrintCode(IndStream &out) {};
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual bool IsTunnel() const {return true;}
  virtual bool IsTunnel(TunType type) const;
  virtual unsigned int NumOutputs() const;
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Tunnel";}
  virtual const DataTypeInfo& DataType(ConnNum num) const;
#if TWOD
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const;
  virtual const SizeList* LocalN(ConnNum num) const;
#endif
#elif DOTENSORS
  virtual const Dim NumDims(ConnNum num) const;
  virtual const SizeList* Len(ConnNum num, Dim dim) const;
  virtual const SizeList* LocalLen(ConnNum num, Dim dim) const;
#endif
  virtual Name GetName(ConnNum num) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual bool Overwrites(const Node *input, ConnNum num) const;
  virtual bool IsSplit() const {return false;}
  virtual bool IsCombine() const {return false;}
  Tunnel* GetRealTunnel();
  const Tunnel* GetRealTunnel() const;
  virtual void MigrateFromOldTun(Tunnel *tun) {}
};

typedef vector<Tunnel*> TunVec;
typedef TunVec::iterator TunVecIter;
typedef TunVec::const_iterator TunVecConstIter;



string TunTypeToStr(TunType type);
