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

#include "loop.h"
#include "loopTunnel.h"

class Combine;

//LoopTunnel for spliting/indxing into a matrix
class Split : public LoopTunnel
{
 public:
  PartDir m_dir;
  bool m_isControlTun;
  bool m_addDir;
  Split();
  Split(PartDir dir, PossTunType type, bool isControl = false);
  ~Split();
  static Node* BlankInst() { return  new Split(LASTPARTDIR,LASTTUNNEL,false);}
  virtual Node* GetNewInst() {return BlankInst(); }
  virtual PossTunnel* GetSetTunnel();
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);
  virtual NodeType GetType() const;
  virtual void SanityCheck();
  virtual unsigned int NumOutputs() const {return GetNumElems(m_dir)+1;}
  virtual bool QuadInUse(Quad quad, bool atEnd) const;
  virtual ClassType GetNodeClass() const {return GetClass();}
  static ClassType GetClass() {return "split";}
  virtual const Sizes* GetM(unsigned int num) const;
  virtual const Sizes* GetN(unsigned int num) const;
  virtual const Sizes* LocalM(unsigned int num) const;
  virtual const Sizes* LocalN(unsigned int num) const;
  void GetSizes(unsigned int num, unsigned int numIters,
		   Size bs, 
		   Size m, Size n,
		   Sizes &ms, Sizes &ns);
  virtual Name GetName(unsigned int num) const;
  virtual Name GetName(unsigned int num, LoopType type) const;
  virtual void PrintVarDeclarations(IndStream &out) const;
  Combine* CreateMatchingCombine(int numArgs, ...);
  bool ValidIter() const;
  unsigned int NumIters(Size bs, Size m, Size n) const;
  unsigned int NumIters(unsigned int iterNum) const;
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);
  unsigned int NumberOfLoopExecs() const;
  void SetAddDir() {m_addDir = true;}
  virtual void StartFillingSizes();
  virtual void ClearSizeCache();
  virtual void AppendSizes(unsigned int execNum, unsigned int numIters);
};
