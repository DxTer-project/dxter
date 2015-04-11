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

#include "realLoop.h"
#include "shadowLoop.h"
#include "loopTunnel.h"
#include "splitBase.h"
#include "LLDLA.h"


class CombineSingleIter;

//LoopTunnel for spliting/indxing into a matrix
class SplitSingleIter : public SplitBase
{
 public:
  bool m_addDir;
  DataTypeInfo m_info;


#if TWOD
  SplitSingleIter();
  SplitSingleIter(PartDir dir, TunType type, bool isControl = false);
#else
  SplitSingleIter(unsigned int partDim, TunType type, bool isControl = false);
#endif

  static Node* BlankInst();
  virtual Node* GetNewInst() {return BlankInst(); }
  virtual Tunnel* GetSetTunnel();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual unsigned int NumOutputs() const;
  virtual bool QuadInUse(Quad quad, bool atEnd) const;
  virtual bool PartInUse(unsigned int partNum) const;
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "split";}
#if TWOD
  virtual const SizeList* GetM(ConnNum num) const;
  virtual const SizeList* GetN(ConnNum num) const;
#if DODM
  virtual const SizeList* LocalM(ConnNum num) const;
  virtual const SizeList* LocalN(ConnNum num) const;
#endif
  void GetSizes(ConnNum num,
		const SizeList *control, Size bs, 
		const SizeList *m, const SizeList *n,
		const SizeList **ms, const SizeList **ns);
#else
  virtual const Dim NumDims(ConnNum num) const;
  virtual const SizeList* Len(ConnNum num, Dim dim) const;
  virtual const SizeList* LocalLen(ConnNum num, Dim dim) const;
#endif
  virtual Name GetName(ConnNum num) const;
  virtual Name GetName(ConnNum num, LoopType type) const;
  virtual void PrintVarDeclarations(BSSize bs, IndStream &out) const;
  CombineSingleIter* CreateMatchingCombine(int numArgs, ...);
  bool ValidIter() const;
#if TWOD
  virtual unsigned int NumIters(Size bs, Size m, Size n) const;
#else
  virtual unsigned int NumIters(Size bs, Size size) const;
#endif
  virtual unsigned int NumIters(unsigned int iterNum) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  virtual unsigned int NumberOfLoopExecs() const;
  void SetAddDir() {m_addDir = true;}

#if DODM
  virtual void UpdateLocalSizes();
#endif
#if DOLLDLA
  virtual string LoopBound();
#endif

  virtual void AddVariables(VarSet &set) const;

  virtual void PrintIncrementAtEndOfLoop(BSSize bs, IndStream &out) const;

#if DOLLDLA
  virtual void BuildDataTypeCache();
  virtual const DataTypeInfo& DataType(ConnNum num) const;
#endif

#if DOLLDLA
  virtual void MigrateFromOldTun(Tunnel *tun);
#endif

  virtual LoopTunnel* GetMatchingOutTun() const;

  string LoopLevel() const;

  virtual void BuildSizes(const SizeList *controlSizes, int stride);

  virtual const SizeList* GetControlSizes() const;
};

