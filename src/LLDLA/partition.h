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

#include "DLANode.h"
#include "node.h"
#include "layers.h"

#if DOLLDLA

class Partition : public DLANode
{
  
 private:
  Dir m_partType;
  Sizes* m_startSizes;
  Sizes* m_endSizes;
  Name startName;
  Name endName;

 public:
  Layer m_layer;
  
  Partition(Layer layer, Dir partType, Size partStart, Size totalSize);

  virtual ~Partition() {}
  inline void SetLayer(Layer layer) {m_layer = layer;}
  inline Layer GetLayer() const {return m_layer;}
  virtual bool IsReadOnly() const {return false;}
  virtual bool CanTrans() const {return false;}

  virtual Name GetName(ConnNum num) const;
  virtual void AddVariables(VarSet &set) const;

  virtual void Prop();

  virtual const Sizes* GetM(ConnNum num) const;
  virtual const Sizes* GetN(ConnNum num) const;

};

#endif //DOLLDLA
