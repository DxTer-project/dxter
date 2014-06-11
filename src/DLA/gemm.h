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

#include "layers.h"
#if (DOBLIS||DOELEM||DOLLDLA)

#include "DLAOp.h"
#include "transform.h"
#include "elemRedist.h"
#include "lowerLayer.h"


#if DOELEM
bool IsDMGemm(const Node *node);
#endif


class Gemm : public DLAOp<3,1>
{
 public:
#if !DOLLDLA
  Trans m_transA, m_transB;
#endif
  Coef m_alpha, m_beta;
  Type m_type;
#if DOBLIS
  Comm m_comm;
#endif

#if !DOLLDLA
  Gemm(Layer layer, Trans transA, Trans transB, Coef alpha, Coef beta, Type type);
#else
  Gemm(Layer layer, Coef alpha, Coef beta, Type type);
#endif

  static Node* BlankInst();
  virtual Node* GetNewInst() { return BlankInst(); }
  virtual void Duplicate(const Node *orig, bool shallow, bool possMerging);

  static ClassType GetClass() {return "Gemm";}
  virtual ClassType GetNodeClass() const {return GetClass();}

  //Ignore these for now
  virtual void FlattenCore(ofstream &out) const;
  virtual void UnflattenCore(ifstream &in, SaveInfo &info);

  virtual NodeType GetType() const;
#if DOELEM
  virtual const DistType& GetDistType(unsigned int num) const;
  virtual bool CanTransposeInputs() const;
  virtual bool DoNotCullDP() const;
#endif
  virtual Phase MaxPhase() const;
  virtual void Prop();
  virtual void PrintCode(IndStream &out);
  
  static Cost GetCost(Layer layer, const Sizes *localDim1, const Sizes *localDim2, const Sizes *localDim3);
#if DOBLIS
  virtual void UpdateInnerPackingMultiple(PackSize size);
  virtual bool IsBLISParallelizable() const;
  virtual void Parallelize(Comm comm);
  virtual bool IsParallel() const;
  virtual bool RemoveParallelization();
  virtual Comm ParallelComm() const {return m_comm;}
#endif
};

#endif
