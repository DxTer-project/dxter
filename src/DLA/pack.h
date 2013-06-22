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

#include "DLANode.h"
#include "DLAOp.h"

enum PackType { PACKROWPANS,
		PACKCOLPANS,
		BADPACK };

enum PackMat { PACKABLOCK,
	       PACKBPANEL,
	       BADPACKMAT };


string PackTypeToStr(PackType type);
string PackSizeToStr(PackSize size);
Size PackSizeToSize(PackSize size);

class Pack : public DLAOp<2,1>
{
 public:
  PackType m_pack;
  unsigned int m_var;
  bool m_scaleAlpha, m_densify, m_invertDiag,
    m_revUpper, m_revLower;
  Comm m_comm;
 Pack(PackType pack, unsigned int var,
      bool scaleAlpha, bool densify, bool invertDiag,
      bool revUpper, bool revLower);
 static Node* BlankInst() {return new Pack(PACKROWPANS, 0, false, false, false, false, false);}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "Pack";}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual bool IsReadOnly() const {return false;}
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual NodeType GetType() const;
  virtual DistType GetDistType(unsigned int num) const;
  virtual Phase MaxPhase() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual Name GetName(unsigned int num) const;
  virtual void SanityCheck();
  virtual unsigned int NumOutputs() const {return 1;}
  inline void Parallelize(Comm comm) {m_comm=comm;}
};

class PackBuff : public DLAOp<1,1>
{
 public:
  PackType m_pack;
  PackSize m_m, m_n;
  PackMat m_packMat;
  Name m_name;
  bool m_densify, m_invertDiag,
    m_revUpper, m_revLower;
  Tri m_tri;
  TriStruct m_triStruct;
  Diag m_diag;
  unsigned int m_parFactor;
  PackBuff(string name, PackType pack,
	   PackMat mat,
	   Tri tri, Diag diag, TriStruct triStruct,
	   bool densify, bool invertDiag, bool revUpper, bool revLower,
	   PackSize mSize, PackSize nSize);
  virtual bool IsReadOnly() const {return true;}
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  static Node* BlankInst() {return new PackBuff("", BADPACK, PACKABLOCK, NOTTRI, NOTTRIDIAG, GEN, false, false, false, false, USEMRSIZE, USEMRSIZE);}
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "PackBuff";}
  virtual NodeType GetType() const { return "packBuff"; }
  virtual DistType GetDistType(unsigned int num) const {return InputDistType(0);}
  virtual Phase MaxPhase() const { return NUMPHASES;}
  virtual void PrintCode(IndStream &out);
  virtual void Prop();
  virtual Name GetName(unsigned int num) const;
  virtual void SanityCheck();
  virtual unsigned int NumOutputs() const {return 1;}
  void UpdateChildrensInnerMultiple(PackSize size);
  virtual bool Overwrites(const Node *input, unsigned int num) const {return false;}
  inline void Parallelize(unsigned int parFactor) {m_parFactor=parFactor;}
};

class LoopInvariantPackBuffMotion : public SingleTrans
{
 public:
  virtual string GetType() const {return "Loop Invariant PackBuff Motion";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class LoopInvariantPackMotion : public SingleTrans
{
 public:
  virtual string GetType() const {return "Loop Invariant Pack Motion";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class CombinePacking : public SingleTrans
{
 public:
  virtual string GetType() const {return "Combine packing";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class CombinePackBuff : public SingleTrans
{
 public:
  virtual string GetType() const {return "Combine Pack Buff";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};


class UnifyPackBuffParams : public SingleTrans
{
 public:
  virtual string GetType() const {return "Unify Pack Buff Params";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

class ReuseTrsmPacking : public SingleTrans
{
 public:
  Layer m_layer;
 ReuseTrsmPacking(Layer layer) : m_layer(layer) {}
  virtual string GetType() const {return "Reuse Trsm packing";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
};

#if DOSQOPHASE
class RenamePackBuff : public SingleTrans
{
 public:
  virtual string GetType() const {return "Rename PackBuff";}
  virtual bool CanApply(const Poss *poss, const Node *node) const;
  virtual void Apply(Poss *poss, Node *node) const;
  string GetNewName(const PackBuff *buff) const;
};
#endif //DOSQOPHASE

bool FindOtherPackBuffs(const Poss *poss, PackMat pack, const Node *ignore);
