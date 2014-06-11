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



#include "DLAOp.h"

#if TWOD
template<unsigned int numIn, unsigned int numOut>
const Sizes* DLAOp<numIn, numOut>::GetM(unsigned int num) const
{
  if (num >= numOut)
    throw;
  return GetInputM(numIn - numOut + num);  
}

template<unsigned int numIn, unsigned int numOut>
   const Sizes* DLAOp<numIn, numOut>::GetN(unsigned int num) const
{
  if (num >= numOut)
    throw;
  return GetInputN(numIn - numOut + num);  
}

template<unsigned int numIn, unsigned int numOut>
   const Sizes* DLAOp<numIn, numOut>::LocalM(unsigned int num) const
{
  if (num >= numOut)
    throw;
  return InputLocalM(numIn - numOut + num);  
}

template<unsigned int numIn, unsigned int numOut>
   const Sizes* DLAOp<numIn, numOut>::LocalN(unsigned int num) const
{
  if (num >= numOut)
    throw;
  return InputLocalN(numIn - numOut + num);  
}

#if DOLLDLA
template<unsigned int numIn, unsigned int numOut>
Stride DLAOp<numIn, numOut>::RowStride(unsigned int num) const
{
  if (num >= numOut)
    throw;
  return InputRowStride(numIn - numOut + num);  
}

template<unsigned int numIn, unsigned int numOut>
   Stride DLAOp<numIn, numOut>::ColStride(unsigned int num) const
{
  if (num >= numOut)
    throw;
  return InputColStride(numIn - numOut + num);  
}
#endif //DOLLDLA

#else
template<unsigned int numIn, unsigned int numOut>
const Sizes* DLAOp<numIn, numOut>::Len(unsigned int num, Dim dim) const
{
  if (num >= numOut)
    throw;
  return InputLen(numIn - numOut + num,dim);  
}

template<unsigned int numIn, unsigned int numOut>
   const Dim DLAOp<numIn, numOut>::NumDims(unsigned int num) const
{
  if (num >= numOut)
    throw;
  return InputNumDims(numIn - numOut + num);  
}

template<unsigned int numIn, unsigned int numOut>
const Sizes* DLAOp<numIn, numOut>::LocalLen(unsigned int num, Dim dim) const
{
  if (num >= numOut)
    throw;
  return InputLocalLen(numIn - numOut + num,dim);  
}

#endif

template<unsigned int numIn, unsigned int numOut>
   Name DLAOp<numIn, numOut>::GetName(unsigned int num) const
{
  if (num >= numOut)
    throw;
  return GetInputName(numIn - numOut + num);  
}

template<unsigned int numIn, unsigned int numOut>
   void DLAOp<numIn, numOut>::Prop()
{
  if (!IsValidCost(m_cost)) {
    for(unsigned int i = 0; i < numIn; ++i)
      Input(i)->Prop();

    if (m_inputs.size() != numIn) {
      cout << "Error on " << GetNodeClass() << endl;
      cout << "Has " << m_inputs.size() << " inputs\n";
      throw;
    }
  }
}

template<unsigned int numIn, unsigned int numOut>
   unsigned int DLAOp<numIn, numOut>::NumOutputs() const
{
  return numOut;
}

template<unsigned int numIn, unsigned int numOut>
bool DLAOp<numIn, numOut>::Overwrites(const Node *input, unsigned int num) const
{
  for(unsigned int i = (numIn-numOut); i < numIn; ++i) {
    const NodeConn *conn = m_inputs[i];
    if (conn->m_n == input && conn->m_num == num)
      return true;
  }
  return false;
}

template<unsigned int numIn, unsigned int numOut>
bool DLAOp<numIn, numOut>::KeepsInputVarLive(Node *input, unsigned int numInArg, unsigned int &numOutArg) const
{
  for(unsigned int i = (numIn-numOut); i < numIn; ++i) {
    const NodeConn *conn = m_inputs[i];
    if (conn->m_n == input && conn->m_num == numInArg) {
      numOutArg = i-(numIn-numOut);
      return true;
    }
  }
  return false;
}

template class DLAOp<2,2>;
template class DLAOp<3,1>;
template class DLAOp<3,2>;
template class DLAOp<2,1>;
template class DLAOp<1,1>;
template class DLAOp<5,1>;
template class DLAOp<5,2>;
template class DLAOp<5,3>;
template class DLAOp<7,7>;
