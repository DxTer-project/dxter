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



#include "DLAOp.h"

template<ConnNum numIn, ConnNum numOut>
const DataTypeInfo& DLAOp<numIn, numOut>::DataType(ConnNum num) const
{
  if (num >= numOut)
    throw;
  return InputDataType(numIn - numOut + num);  
}

#if TWOD
template<ConnNum numIn, ConnNum numOut>
const SizeList* DLAOp<numIn, numOut>::GetM(ConnNum num) const
{
  if (num >= numOut)
    throw;
  return GetInputM(numIn - numOut + num);  
}

template<ConnNum numIn, ConnNum numOut>
   const SizeList* DLAOp<numIn, numOut>::GetN(ConnNum num) const
{
  if (num >= numOut)
    throw;
  return GetInputN(numIn - numOut + num);  
}

#if DODM
template<ConnNum numIn, ConnNum numOut>
   const SizeList* DLAOp<numIn, numOut>::LocalM(ConnNum num) const
{
  if (num >= numOut)
    throw;
  return InputLocalM(numIn - numOut + num);  
}

template<ConnNum numIn, ConnNum numOut>
   const SizeList* DLAOp<numIn, numOut>::LocalN(ConnNum num) const
{
  if (num >= numOut)
    throw;
  return InputLocalN(numIn - numOut + num);  
}
#endif //DODM


#else
template<ConnNum numIn, ConnNum numOut>
const SizeList* DLAOp<numIn, numOut>::Len(ConnNum num, Dim dim) const
{
  if (num >= numOut)
    throw;
  return InputLen(numIn - numOut + num,dim);  
}

template<ConnNum numIn, ConnNum numOut>
   const Dim DLAOp<numIn, numOut>::NumDims(ConnNum num) const
{
  if (num >= numOut)
    throw;
  return InputNumDims(numIn - numOut + num);  
}

template<ConnNum numIn, ConnNum numOut>
const SizeList* DLAOp<numIn, numOut>::LocalLen(ConnNum num, Dim dim) const
{
  if (num >= numOut)
    throw;
  return InputLocalLen(numIn - numOut + num,dim);  
}

#endif

template<ConnNum numIn, ConnNum numOut>
   Name DLAOp<numIn, numOut>::GetName(ConnNum num) const
{
  if (num >= numOut)
    throw;
  return GetInputName(numIn - numOut + num);  
}

template<ConnNum numIn, ConnNum numOut>
   void DLAOp<numIn, numOut>::Prop()
{
  if (!IsValidCost(m_cost)) {
    for(ConnNum i = 0; i < numIn; ++i)
      Input(i)->Prop();

    if (m_inputs.size() != numIn) {
      cout << "Error on " << GetNodeClass() << endl;
      cout << "Has " << m_inputs.size() << " inputs\n";
      throw;
    }
  }
}

template<ConnNum numIn, ConnNum numOut>
   unsigned int DLAOp<numIn, numOut>::NumOutputs() const
{
  return numOut;
}

template<ConnNum numIn, ConnNum numOut>
bool DLAOp<numIn, numOut>::Overwrites(const Node *input, ConnNum num) const
{
  for(ConnNum i = (numIn-numOut); i < numIn; ++i) {
    const NodeConn *conn = m_inputs[i];
    if (conn->m_n == input && conn->m_num == num)
      return true;
  }
  return false;
}
/*
template<ConnNum numIn, ConnNum numOut>
bool DLAOp<numIn, numOut>::KeepsInputVarLive(Node *input, ConnNum numInArg, ConnNum &numOutArg) const
{
  for(ConnNum i = (numIn-numOut); i < numIn; ++i) {
    const NodeConn *conn = m_inputs[i];
    if (conn->m_n == input && conn->m_num == numInArg) {
      numOutArg = i-(numIn-numOut);
      return true;
    }
  }
  return false;
}
*/

template class DLAOp<2,2>;
template class DLAOp<3,1>;
template class DLAOp<4,1>;
template class DLAOp<3,2>;
template class DLAOp<2,1>;
template class DLAOp<1,1>;
template class DLAOp<5,1>;
template class DLAOp<5,2>;
template class DLAOp<5,3>;
template class DLAOp<7,7>;
