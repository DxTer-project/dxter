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



#include "elemRedist.h"
#include "math.h"
#include "helperNodes.h"

#if DOELEM

//void GetLocalSize(DistType dist, Size m, Size n, Size &localM, Size &localN);

bool CanTrans(DistType src, DistType dest, bool tryOpposite)
{
  if (src == D_MC_MR && dest == D_STAR_MR)
    return true;
  if (src == D_MR_MC && dest == D_STAR_MC)
    return true;
  if (src == D_VC_STAR && dest == D_MC_STAR)
    return true;
  if (src == D_VR_STAR && dest == D_MR_STAR)
    return true;

  if (tryOpposite)
    return CanTrans(dest, src, false);

  return false;
}

Trans SwitchTrans(Trans trans, Type type)
{
  if (trans == NORMAL) {
    if (type == REAL)
      return TRANS;
    else
      return CONJTRANS;
  }
  else if (trans == TRANS)
    return NORMAL;
  else if (trans == CONJTRANS)
    return NORMAL;
  throw;
}

DistType TransType(DistType dist, Trans trans)
{
  switch(dist) {
  case (D_STAR_MR):
    if (trans == TRANS)
      return D_MR_STAR_T;
    else
      return D_MR_STAR_H;
  case (D_STAR_MC):
    if (trans == TRANS)
      return D_MC_STAR_T;
    else
      return D_MC_STAR_H;
  case (D_MR_STAR):
    if (trans == TRANS)
      return D_STAR_MR_T;
    else
      return D_STAR_MR_H;
  case (D_MC_STAR):
    if (trans == TRANS)
      return D_STAR_MC_T;
    else
      return D_STAR_MC_H;
  default:
    throw;
  }
}


Trans UpdateTrans(Trans trans, DistType dist)
{
  Trans inputTrans = DistTransType(dist);
  switch (trans)
  {
    case (NORMAL):
      return inputTrans;
    case (TRANS):
      if (inputTrans == TRANS)
        return NORMAL;
      else if (inputTrans == CONJTRANS) {
        throw;
      }
      else
        return TRANS;
    case (CONJTRANS):
      if (inputTrans == TRANS) {
        throw;
      }
      else if (inputTrans == CONJTRANS)
        return NORMAL;
      else
        return CONJTRANS;
    default:
      throw;
      return NORMAL;
      break;
  }
}

DistType GetNonTrans(DistType dist)
{
  switch(dist) {
    case D_MC_STAR_T:
    case D_MC_STAR_H:
      return D_STAR_MC;
    case D_STAR_MC_T:
    case D_STAR_MC_H:
      return D_MC_STAR;
    case D_MR_STAR_T:
    case D_MR_STAR_H:
      return D_STAR_MR;
    case D_STAR_MR_T:
    case D_STAR_MR_H:
      return D_MR_STAR;
    default:
      throw;
  }
  
  
}

Trans DistTransType(DistType dist)
{
  switch(dist) {
    case D_MC_STAR_T:
      return TRANS;
    case D_MC_STAR_H:
      return CONJTRANS;
    case D_STAR_MC_T:
      return TRANS;
    case D_STAR_MC_H:
      return CONJTRANS;
    case D_MR_STAR_T:
      return TRANS;
    case D_MR_STAR_H:
      return CONJTRANS;
    case D_STAR_MR_T:
      return TRANS;
    case D_STAR_MR_H:
      return CONJTRANS;
    case D_MR_MC_H:
      return CONJTRANS;
    case D_MR_MC_T:
      return TRANS;
    case D_MC_MR_H:
      return CONJTRANS;
    case D_MC_MR_T:
      return TRANS;
    case D_STAR_STAR:
    case D_MC_MR:
    case D_MR_MC:
    case D_MC_STAR:
    case D_STAR_MC:
    case D_MR_STAR:
    case D_STAR_MR:
    case D_VC_STAR:
    case D_STAR_VC:
    case D_VR_STAR:
    case D_STAR_VR:
      return NORMAL;
  default:
    throw;
  }
}

bool IsTransType(DistType dist)
{
  switch(dist) {
    case D_MC_STAR_T:
    case D_MC_STAR_H:
    case D_STAR_MC_T:
    case D_STAR_MC_H:
    case D_MR_STAR_T:
    case D_MR_STAR_H:
    case D_STAR_MR_T:
    case D_STAR_MR_H:
      return true;
    case D_STAR_STAR:
    case D_MC_MR:
    case D_MR_MC:
    case D_MC_STAR:
    case D_STAR_MC:
    case D_MR_STAR:
    case D_STAR_MR:
    case D_VC_STAR:
    case D_STAR_VC:
    case D_VR_STAR:
    case D_STAR_VR:
      return false;
    default:
      throw;
  }
}

template<>
Cost GetRedistCost<D_VC_STAR>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
    case(D_MC_STAR):
      return CopyCost(ceil((double)m/P), n, C, ONE);
    case(D_MC_MR):
      return AllToAll(C * ceil((double)m/P * (double)n/C), C) + 
      (C*CopyCost(ceil((double)m/P), ceil((double)n/C), C, ONE)) + 
      (C*CopyCost(ceil((double)m/P), ceil((double)n/C), ONE, ONE));
    case(D_VR_STAR):
      return SendRecv(ceil((double)m / P) * n) + CopyCost(ceil((double)m/P), n, ONE, ONE) + CopyCost(ceil((double)m/P), n, ONE, ONE);
  case (D_STAR_MC_T):
  case (D_STAR_MC_H):
    return INFINITY;
  case (D_MR_STAR):
      return CopyCost(ceil((double)m/P), n, R, ONE); + GetRedistCost<D_VC_STAR>(D_VR_STAR,m,n);
    default:
      throw;
  }
}

template<>
Cost GetRedistCost<D_MC_STAR>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
    case(D_VC_STAR):
      return AllGather(C * ceil((double)m / P) * n, C) 
	+ CopyCost(ceil((double)m/P), n, ONE, ONE) 
	+ (C * CopyCost(ceil((double)m/P), n, ONE, C));
    case (D_VR_STAR):
      return GetRedistCost<D_VC_STAR>(D_VR_STAR,m,n) + GetRedistCost<D_MC_STAR>(D_VC_STAR,m,n);
    case(D_MC_MR):
      return AllGather(C * ceil((double)m / R * (double)n / C), C) + 
	CopyCost(ceil((double)m/R), ceil((double)n/C), ONE, ONE) + 
	(C * CopyCost(ceil((double)m/R), ceil((double)n/C), ONE, ONE)); 
    default:
      throw;
  }
}

template<>
Cost GetRedistCost<D_STAR_VR>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
    case (D_MC_MR):
      return AllToAll(R*ceil((double)m/R * (double)n/P),R) 
      + (R * CopyCost(ceil((double)m/R),ceil((double)n/P),ONE,ONE)) 
      + (R * CopyCost(ceil((double)m/R),ceil((double)n/P),ONE,R));
    case (D_STAR_VC):
      return SendRecv(m*ceil((double)n/P)) + CopyCost(m,ceil((double)n/P),ONE,ONE) + CopyCost(m,ceil((double)n/P),ONE,ONE);
    case (D_STAR_MR):
      return CopyCost(m,ceil((double)n/P), ONE, ONE, true, false);
    case (D_MR_MC):
      //Inlining code to remove a cicular dependency on explicitly instantiated templates
      //    return GetRedistCost<D_STAR_VC>(D_MR_MC,m,n) + GetRedistCost<D_STAR_VR>(D_STAR_VC,m,n);
      return GetRedistCost<D_STAR_VR>(D_MC_MR,m,n) + SendRecv(m * ceil((double)n/P)) + CopyCost(m,ceil((double)n/P),ONE,ONE) + CopyCost(m,ceil((double)n/P),ONE,ONE)
      + GetRedistCost<D_STAR_VR>(D_STAR_VC,m,n);
    case (D_MR_STAR_T):
      return CopyCost(m, ceil((double)n/P), ceil((double)n/C), ONE, true, false);
    case (D_MR_STAR_H):
      return CopyCost(m, ceil((double)n/P), ceil((double)n/C), ONE, true, false) + GAMMA * m * ceil((double)n/P);
    default: 
      throw;
  }
}

template<>
Cost GetRedistCost<D_STAR_VC>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
    case(D_STAR_MC):
      return CopyCost(m, ceil((double)n/P), ONE, ONE);
    case (D_STAR_VR):
      return SendRecv(m * ceil((double)n/P)) + CopyCost(m,ceil((double)n/P),ONE,ONE) + CopyCost(m,ceil((double)n/P),ONE,ONE);
    case (D_MC_MR):
      return GetRedistCost<D_STAR_VR>(D_MC_MR,m,n) + GetRedistCost<D_STAR_VC>(D_STAR_VR,m,n);
    case (D_STAR_STAR):
      return CopyCost(m, ceil((double)n/P), ONE, ONE);
    default:
      throw;
  }
}

template<>
Cost GetRedistCost<D_VR_STAR>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
    case (D_MR_STAR):
      return CopyCost(ceil((double)m/P), n, R, ONE);
    case (D_MC_MR):
      return GetRedistCost<D_VC_STAR>(D_MC_MR,m,n) + GetRedistCost<D_VR_STAR>(D_VC_STAR,m,n);
    case(D_VC_STAR):
      return SendRecv(ceil((double)m / P) * n) + CopyCost(ceil((double)m/P), n, ONE, ONE) + CopyCost(ceil((double)m/P), n, ONE, ONE);
    case (D_MR_MC):
      return (R * CopyCost(ceil((double)m/P),ceil((double)n/R),R,ONE)) 
      + AllToAll(R*ceil((double)m/P * (double)n/R),R) + (R* CopyCost(ceil((double)m/P),ceil((double)n/R),ONE,ONE));
    default:
      throw;
  }
}

template<>
Cost GetRedistCost<D_VC_STAR_T>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case (D_MR_STAR_T):
    return GetRedistCost<D_VR_STAR>(D_MR_STAR, n, m) 
      + GetRedistCost<D_VC_STAR>(D_VR_STAR, n, m);
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_VC_STAR_H>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case (D_MR_STAR_H):
    return GetRedistCost<D_VC_STAR_T>(D_MR_STAR_T, m, n);
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_MC_MR>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
    case(D_VC_STAR):
      return AllToAll(C*ceil((double)m/P)*ceil((double)n/C),C)
	+ (C * CopyCost(ceil((double)m/P), ceil((double)n/C), ONE, ONE)) 
	+ (C * CopyCost(ceil((double)m/P), ceil((double)n/C), ONE, C));
    case(D_VR_STAR):
      return GetRedistCost<D_VC_STAR>(D_VR_STAR,m,n) + GetRedistCost<D_MC_MR>(D_VC_STAR,m,n);
    case(D_MC_STAR):
      return CopyCost(ceil((double)m/R), ceil((double)n/C), ONE, ONE);
    case (D_MR_STAR):
      return GetRedistCost<D_VR_STAR>(D_MR_STAR,m,n) + GetRedistCost<D_VC_STAR>(D_VR_STAR,m,n) + GetRedistCost<D_MC_MR>(D_VC_STAR,m,n);
    case(D_STAR_STAR):
      return CopyCost(ceil((double)m/R), ceil((double)n/C), R, ONE);
    case(D_STAR_MC):
      return GetRedistCost<D_STAR_VC>(D_STAR_MC,m,n) + 
	GetRedistCost<D_STAR_VR>(D_STAR_VC,m,n) + GetRedistCost<D_MC_MR>(D_STAR_VR,m,n);
    case (D_STAR_MC_T):
    case (D_STAR_MC_H):
      return CopyCost(ceil((double)m/R), ceil((double)n/C), n, ONE);
    case (D_STAR_VC):
      return GetRedistCost<D_STAR_VR>(D_STAR_VC,m,n) + GetRedistCost<D_MC_MR>(D_STAR_VR,m,n);
    case (D_STAR_VR):
      return AllToAll(R*ceil((double)m/R*(double)n/P),R) 
      + (R * CopyCost(ceil((double)m/R),ceil((double)n/P),R,ONE)) 
      + (R * CopyCost(ceil((double)m/R),ceil((double)n/P),ONE,ONE));
    case (D_STAR_MR):
      return CopyCost(ceil((double)m/R),ceil((double)n/C),R,ONE);
    case (D_MR_MC): 
      return GetRedistCost<D_VR_STAR>(D_MR_MC, m, n) + GetRedistCost<D_VC_STAR>(D_VR_STAR,m,n) + GetRedistCost<D_MC_MR>(D_VC_STAR,m,n); 
    case (D_MR_STAR_T):
      return CopyCost(ceil((double)m/R),ceil((double)n/C),R,ONE);
    case (D_MR_STAR_H):
      return CopyCost(ceil((double)m/R),ceil((double)n/C),R,ONE) + GAMMA * ceil((double)m/R * (double)n/C);
  case (D_MR_MC_H):
  case (D_MR_MC_T):
    return CopyCost(ceil((double)n/C),ceil((double)m/R),ONE,ceil((double)m/R));
  default:
    throw;
  }
}


template<>
Cost GetRedistCost<D_MR_MC>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_VR_STAR):
    return AllToAll(R*ceil((double)m/P*(double)n/R),R) + R*(CopyCost(ceil((double)m/P),ceil((double)n/R),ONE,ONE)) + R*CopyCost(ceil((double)m/P),ceil((double)n/R),ONE,R);
  case(D_MC_MR):
    return GetRedistCost<D_VC_STAR>(srcType,m,n) + GetRedistCost<D_VR_STAR>(D_VC_STAR,m,n) + GetRedistCost<D_MR_MC>(D_VR_STAR,m,n);
  case (D_STAR_VC):
    return AllToAll(C*ceil((double)m/C*(double)n/P),C) + (C*CopyCost(ceil((double)m/C),ceil((double)n/P),C,ONE)) + (C*CopyCost(ceil((double)m/C),ceil((double)n/P),ONE,ONE));
    break;
  default:
    throw;
  }
}


template<>
Cost GetRedistCost<D_MR_STAR>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_VR_STAR):
    return AllGather(R * ceil((double)m / P) * n, R) + CopyCost(n, ONE, ONE, ONE) + (R * CopyCost(ceil((double)m/P), n, ceil((double)m/P), R));
  case (D_VC_STAR):
    return GetRedistCost<D_VR_STAR>(D_VC_STAR,m,n)
      + GetRedistCost<D_MR_STAR>(D_VR_STAR,m,n);
  case(D_MR_MC):
    return AllGather(R * ceil((double)m / C * (double)n / R), R) + CopyCost(ceil((double)m/C), ceil((double)n/R), ONE, ONE) + (R * CopyCost(ceil((double)n/R), R, ONE, ONE));
  case (D_MC_MR):
    return GetRedistCost<D_VC_STAR>(D_MC_MR, m, n)
      + GetRedistCost<D_VR_STAR>(D_VC_STAR, m, n)
      + GetRedistCost<D_MR_STAR>(D_VR_STAR, m, n);
  case (D_STAR_STAR):
    return CopyCost(ceil((double)m/C),n,C,ONE);
  case (D_STAR_VR):
    return GetRedistCost<D_STAR_VC>(D_STAR_VR,m,n) + GetRedistCost<D_MR_MC>(D_STAR_VC,m,n) + GetRedistCost<D_STAR_VR>(D_MR_MC,m,n);
  case (D_STAR_VC):
    return GetRedistCost<D_MR_MC>(D_STAR_VC,m,n) + GetRedistCost<D_MR_STAR>(D_MR_MC,m,n);
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_MR_STAR_T>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_MC_MR):
    return AllGather(R*ceil((double)n/C * (double)m/R),R) 
      + CopyCost(ceil((double)n/C),ceil((double)m/R),ceil((double)m/R),ONE, true, false) 
      + (R * CopyCost(ceil((double)n/C),ceil((double)m/R),ONE,ONE, false, true));
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_MR_STAR_H>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_MC_MR):
    return GetRedistCost<D_MR_STAR_T>(srcType, m, n) + GAMMA * ceil((double)n/C * (double)m/R);
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_MR_MC_T>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_MC_MR_T):
    return GetRedistCost<D_MR_MC>(D_MC_MR, n, m);
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_MR_MC_H>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_MC_MR_H):
    return GetRedistCost<D_MR_MC_T>(D_MC_MR_T, m, n);
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_STAR_STAR>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_MC_MR):
    return AllGather(m*n, P) + CopyCost(ceil((double)m/R), ceil((double)m/C), ONE, ONE) + (C * CopyCost(ceil((double)m/R), ceil((double)n/C), ONE, R));
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_STAR_MC_T>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_VC_STAR):
    return AllGather(C * n * ceil((double)m / P), C) 
      + CopyCost(n,ceil((double)m/P),ceil((double)m/P),ONE) 
      + (C*CopyCost(n,ceil((double)m/P),ONE,ONE));
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_STAR_MC_H>(DistType srcType, Size m, Size n)
{
  return GetRedistCost<D_STAR_MC_T>(srcType, m, n) + GAMMA * n * ceil((double)m / C);
  //BAM + conj cost;
}

template<>
Cost GetRedistCost<D_STAR_MR>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case (D_STAR_VR):
    return AllGather(R*m*ceil((double)n/P),R) 
      + CopyCost(m,ceil((double)n/P),ONE,ONE) 
      + (R * CopyCost(m,ceil((double)n/P),ONE,ONE));
  case(D_MC_MR):
    return AllGather(R*ceil((double)m/R*(double)n/C),R) 
      + CopyCost(ceil((double)m/R),ceil((double)n/C),ONE,ONE, false, true) 
      + (R * CopyCost(ceil((double)m/R),ceil((double)n/C),ONE,R, false, true));
  case(D_STAR_STAR):
    return CopyCost(m, ceil((double)n/C), ONE, ONE);
  case (D_STAR_VC):
    return GetRedistCost<D_STAR_VR>(D_STAR_VC,m,n) + GetRedistCost<D_STAR_MR>(D_STAR_VR,m,n);
  default:
    throw;
  }
}


template<>
Cost GetRedistCost<D_STAR_MR_T>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_VR_STAR):
    return AllGather(R * n * ceil((double)m / P), R) + CopyCost(n,ceil((double)m/P),ceil((double)m/P),ONE) + (R*CopyCost(n,ceil((double)m/P),ONE,ONE));
  default:
    throw;
  }
}

template<>
Cost GetRedistCost<D_STAR_MR_H>(DistType srcType, Size m, Size n)
{
  return GetRedistCost<D_STAR_MR_T>(srcType, m, n) + GAMMA * n * ceil((double)m/P);
  //BAM + conj cost;
}


template<>
Cost GetRedistCost<D_STAR_MC>(DistType srcType, Size m, Size n)
{
  switch(srcType) {
  case(D_MC_MR):
    return GetRedistCost<D_STAR_VR>(srcType,m,n) 
      + GetRedistCost<D_STAR_VC>(D_STAR_VR,m,n) 
      + GetRedistCost<D_STAR_MC>(D_STAR_VC,m,n);
  case (D_STAR_VC):
    return AllGather(C*m*ceil((double)n/P),C) + CopyCost(m,ceil((double)n/P),ONE,ONE) + (C * CopyCost(m,ceil((double)n/P),ONE,ONE));
  case (D_STAR_VR):
    return GetRedistCost<D_STAR_VC>(D_STAR_VR,m,n) + GetRedistCost<D_STAR_MC>(D_STAR_VC,m,n);
  case (D_VR_STAR):
    return GetRedistCost<D_MR_MC>(D_VR_STAR,m,n) + GetRedistCost<D_STAR_MC>(D_MR_MC,m,n);
  case (D_VC_STAR):
    return GetRedistCost<D_VR_STAR>(D_VC_STAR,m,n) + GetRedistCost<D_MR_MC>(D_VR_STAR,m,n) + GetRedistCost<D_STAR_MC>(D_MR_MC,m,n);
  case (D_MR_MC):
    return AllGather(C * ceil((double)m/C*(double)n/R),C) + CopyCost(ceil((double)m/C), ceil((double)n/R), ONE, ONE) + (C*CopyCost(ceil((double)m/C),ceil((double)n/R),ONE,C));
  case (D_STAR_MR):
    return GetRedistCost<D_STAR_VR>(D_STAR_MR, m, n) + GetRedistCost<D_STAR_VC>(D_STAR_VR, m, n) + GetRedistCost<D_STAR_MC>(D_STAR_VC, m, n);
  case (D_STAR_STAR):
    return CopyCost(m, ceil((double)n/R), ONE, ONE);
  case (D_VC_STAR_H):
  case (D_VC_STAR_T):
    return CopyCost(m, ceil((double)n/P), ceil((double)n/P), ONE) 
      + AllGather(C * m * ceil((double)n/P), C)
      + C * (CopyCost(m, ceil((double)n/P), ONE, ONE));
  default:
    throw;
  }
}

//Consider the distribution of the ORIGINAL data,
// not of the transposed data
OneDDistType GetColDist(DistType dist)
{
  switch (dist)
    {
    case (D_STAR_STAR):
    case (D_STAR_MC):
    case (D_MC_STAR_T):
    case (D_MC_STAR_H):
    case (D_STAR_MR):
    case (D_MR_STAR_T):
    case (D_MR_STAR_H):
    case (D_STAR_VR):
    case (D_STAR_VC):
      return STAR;
    case (D_MC_MR):
    case (D_MC_STAR):
    case (D_STAR_MC_T):
    case (D_STAR_MC_H):
      return MC;
    case (D_MR_MC):
    case (D_MR_STAR):
    case (D_STAR_MR_T):
    case (D_STAR_MR_H):
      return MR;
    case (D_VC_STAR):
      return VC;
    case (D_VR_STAR):
      return VR;
    case (UNKNOWN):
    case (D_LASTDIST):
    default:
      throw;
    }
}

//Consider the distribution of the ORIGINAL data,
// not of the transposed data
OneDDistType GetRowDist( DistType dist)
{
  switch (dist)
    {
    case (D_STAR_STAR):
    case (D_MC_STAR):
    case (D_STAR_MC_T):
    case (D_STAR_MC_H):
    case (D_MR_STAR):
    case (D_STAR_MR_T):
    case (D_STAR_MR_H):
    case (D_VC_STAR):
    case (D_VR_STAR):
      return STAR;
    case (D_MC_MR):
    case (D_STAR_MR):
    case (D_MR_STAR_T):
    case (D_MR_STAR_H):
      return MR;
    case (D_MC_STAR_T):
    case (D_MC_STAR_H):
    case (D_MR_MC):
    case (D_STAR_MC):
      return MC;
    case (D_STAR_VC):
      return VC;
    case (D_STAR_VR):
      return VR;
    case (UNKNOWN):
    case (D_LASTDIST):
    default:
      throw;
    }
}

RemoveWastedRedist::RemoveWastedRedist(DistType destType)
  : m_destType(destType)
{
  m_type = "Removing wasted redist from something to " + 
    DistTypeToStr(destType);
}


bool RemoveWastedRedist::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    return false;
  RedistNode *redistNode = (RedistNode*)node;
  if (redistNode->GetDistType(0) != m_destType)
    return false;
  if (node->m_children.size() == 0)
    throw;
  while(redistNode->Input(0) 
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
    {
      redistNode = (RedistNode*)redistNode->Input(0);
      if (redistNode->GetDistType(0) == m_destType)
	return true;
      for(unsigned int i = 0; i < redistNode->m_children.size(); ++i) {
	Node *tmp = redistNode->Child(i);
	if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
	  if (((RedistNode*)tmp)->m_destType == m_destType)
	    return true;
	}
      }
    }
  if (redistNode->Input(0)
      && (((DLANode*)(redistNode->Input(0)))->GetDistType(redistNode->InputConnNum(0)) == m_destType))
    return true;
  return false;
}

void RemoveWastedRedist::Apply(Node *node) const
{
  RedistNode *redistNode = (RedistNode*)node;
  while(redistNode->Input(0) 
        && (redistNode->Input(0)->GetNodeClass() == RedistNode::GetClass()))
    {
      redistNode = (RedistNode*)redistNode->Input(0);
      if (redistNode->GetDistType(0) == m_destType) {
	node->RedirectChildren(redistNode, 0);
	node->m_poss->DeleteChildAndCleanUp(node);
	return;
      }
      for(unsigned int i = 0; i < redistNode->m_children.size(); ++i) {
	Node *tmp = redistNode->Child(i);
	if (tmp != node && tmp->GetNodeClass() == RedistNode::GetClass()) {
	  if (((RedistNode*)tmp)->m_destType == m_destType) {
	    node->RedirectChildren(tmp, 0);
	    node->m_poss->DeleteChildAndCleanUp(node);
	    return;
	  }
	}
      }
    }
  if (redistNode->Input(0)
      && (((DLANode*)redistNode->Input(0))->GetDistType(redistNode->InputConnNum(0)) == m_destType))
    {
      node->RedirectChildren(redistNode->Input(0), redistNode->InputConnNum(0));
      node->m_poss->DeleteChildAndCleanUp(node);
      return;
    }
  throw;
}

template<DistType SrcType, DistType DestType>
ExpandRedistribution<SrcType,DestType>::ExpandRedistribution()
{
  m_type = "Expanding redistribution " + DistTypeToStr(SrcType) + " -> " +
    DistTypeToStr(DestType);
}

template<DistType SrcType, DistType DestType>
bool ExpandRedistribution<SrcType, DestType>::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    return false;
  RedistNode *redistNode = (RedistNode*)node;
  if (redistNode->GetDistType(0) != DestType)
    return false;
  Node *input = redistNode->Input(0);
  if (!input 
      || (((DLANode*)input)->GetDistType(redistNode->InputConnNum(0)) != SrcType))
    return false;
  return true;
}


template<DistType SrcType, DistType DestType>
bool ExpandRedistribution<SrcType, DestType>::WorthApplying(const Node *node) const
{
  if (CanUseUp(node)) {
    return true;
  }
  if (ExpansionHasPossibleTrans()) {
    for (unsigned int i = 0; i < node->m_children.size(); ++i) {
      const DLANode *child = (DLANode*)(node->Child(i));
      if (child->CanTransposeInputs()) {
        return true;      
      }
    }
    return false;
  }
  return false;
}


template<>
void ExpandRedistribution<D_MC_MR,D_VR_STAR>::Apply(Node *node) const
{
  RedistNode *node1 = new RedistNode(D_VC_STAR);
  RedistNode *node2 = new RedistNode(D_VR_STAR);
  node2->AddInput(node1,0);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node->m_poss->AddNodes(2,node1,node2);
  node->RedirectChildren(node2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

template<>
bool ExpandRedistribution<D_MC_MR,D_VR_STAR>::ExpansionHasPossibleTrans() const
{
    //  return CanTrans(D_MC_MR, D_VC_STAR) || CanTrans(D_VC_STAR, D_MC_STAR);
  //  return CanTrans(D_VC_STAR, D_MC_STAR);
  return false;
}


template<>
void ExpandRedistribution<D_MC_MR,D_MC_STAR>::Apply(Node *node) const
{
  RedistNode *node1 = new RedistNode(D_VC_STAR);
  RedistNode *node2 = new RedistNode(D_MC_STAR);
  node2->AddInput(node1,0);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node->m_poss->AddNodes(2,node1,node2);
  node->RedirectChildren(node2, 0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

template<>
bool ExpandRedistribution<D_MC_MR,D_MC_STAR>::ExpansionHasPossibleTrans() const
{
  //  return CanTrans(D_MC_MR, D_VC_STAR) || CanTrans(D_VC_STAR, D_MC_STAR);
  return CanTrans(D_VC_STAR, D_MC_STAR);
}

template<>
void ExpandRedistribution<D_MC_MR,D_STAR_MR>::Apply(Node *node) const
{
  RedistNode *node1 = new RedistNode(D_STAR_VR);
  RedistNode *node2 = new RedistNode(D_STAR_MR);
  node2->AddInput(node1,0);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node->m_poss->AddNodes(2,node1,node2);
  node->RedirectChildren(node2,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

template<>
bool ExpandRedistribution<D_MC_MR,D_STAR_MR>::ExpansionHasPossibleTrans() const
{
  //  return CanTrans(D_MC_MR, D_STAR_VR) || CanTrans(D_STAR_VR,D_STAR_MR);
  return CanTrans(D_STAR_VR,D_STAR_MR);
}

template<>
void ExpandRedistribution<D_MC_MR,D_MR_STAR>::Apply(Node *node) const
{
  RedistNode *node1 = new RedistNode(D_VC_STAR);
  RedistNode *node2 = new RedistNode(D_VR_STAR);
  RedistNode *node3 = new RedistNode(D_MR_STAR);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node1,0);
  node3->AddInput(node2,0);
  node->m_poss->AddNodes(3,node1,node2,node3);
  node->RedirectChildren(node3,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

template<>
bool ExpandRedistribution<D_MC_MR,D_MR_STAR>::ExpansionHasPossibleTrans() const
{
  //  return CanTrans(D_MC_MR, D_VC_STAR) || CanTrans(D_VC_STAR, D_VR_STAR) 
  //    || CanTrans(D_VR_STAR, D_MR_STAR);
  return CanTrans(D_VR_STAR, D_MR_STAR);
}

template<>
void ExpandRedistribution<D_MC_MR,D_STAR_MC>::Apply(Node *node) const
{
  RedistNode *node1 = new RedistNode(D_STAR_VR);
  RedistNode *node2 = new RedistNode(D_STAR_VC);
  RedistNode *node3 = new RedistNode(D_STAR_MC);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node2->AddInput(node1,0);
  node3->AddInput(node2,0);
  node->m_poss->AddNodes(3,node1,node2,node3);
  node->RedirectChildren(node3,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

template<>
bool ExpandRedistribution<D_MC_MR,D_STAR_MC>::ExpansionHasPossibleTrans() const
{
  //  return CanTrans(D_MC_MR, D_STAR_VR) || CanTrans(D_STAR_VR, D_STAR_VC) 
  //    || CanTrans(D_STAR_VC, D_STAR_MC);
  return CanTrans(D_STAR_VC, D_STAR_MC);
}

template<>
void ExpandRedistribution<D_VC_STAR,D_MR_STAR>::Apply(Node *node) const
{
  RedistNode *node1 = new RedistNode(D_VR_STAR);
  RedistNode *node2 = new RedistNode(D_MR_STAR);
  node2->AddInput(node1,0);
  node1->AddInput(node->Input(0),node->InputConnNum(0));
  node->m_poss->AddNodes(2,node1,node2);
  node->RedirectChildren(node2,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}


template<>
bool ExpandRedistribution<D_VC_STAR,D_MR_STAR>::ExpansionHasPossibleTrans() const
{
  //  return CanTrans(D_VC_STAR,D_VR_STAR) || CanTrans(D_VR_STAR,D_MR_STAR);
  return CanTrans(D_VR_STAR,D_MR_STAR);
}


CombineRedistribs::CombineRedistribs(DistType srcType, DistType destType) 
  : m_srcType(srcType), m_destType(destType)
{
  m_type = "Combine redistributions " + DistTypeToStr(m_srcType) + " -> " +
    DistTypeToStr(m_destType);
}

bool CombineRedistribs::CanApply(const Node *node) const
{
  const static ClassType classType = RedistNode::GetClass();
  if (((RedistNode*)node)->m_destType != m_destType)
    return false;
  const Node *parent = node->Input(0);
  if (!parent
      || (((DLANode*)parent)->GetDistType(node->InputConnNum(0)) != m_srcType))
    return false;
  NodeConnVecConstIter iter = parent->m_children.begin();
  for( ; iter != parent->m_children.end(); ++iter) {
    DLANode *output = (DLANode*)((*iter)->m_n);
    if (output != node
        && (output->GetNodeClass() == classType)
        && (((RedistNode*)output)->m_destType == m_destType)
        && (output->InputConnNum(0) == node->InputConnNum(0))) {
      return true;
    }
  }
  return false;
}

void CombineRedistribs::Apply(Node *node) const
{
  Node *parent = node->Input(0);
  NodeConnVecIter iter = parent->m_children.begin();
  for( ; iter != parent->m_children.end(); ++iter) {
    DLANode *output = (DLANode*)(*iter)->m_n;
    if (output != node
        && (output->GetNodeClass() == RedistNode::GetClass())
        && (output->GetDistType(0) == m_destType)
        && (output->InputConnNum(0) == node->InputConnNum(0))) {
      output->RedirectChildren(node,0);
      output->m_poss->DeleteChildAndCleanUp(output);
      return;
    }
  }
}

RemoveNOPRedistribs::RemoveNOPRedistribs(DistType type)
  : m_distribType(type)
{
  m_type = "Remove NOP Redistribution " + DistTypeToStr(m_distribType);
}

bool RemoveNOPRedistribs::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    return false;
  DLANode *ddla = (DLANode*)node;
  if (ddla->GetDistType(0) != m_distribType) {
    return false;
  }
  Node *parent = ddla->Input(0);
  if (!parent
      || (((DLANode*)parent)->GetDistType(node->InputConnNum(0)) != m_distribType)) {
    return false;
  }
  return true;
}

void RemoveNOPRedistribs::Apply(Node *node) const
{
  Node *parent = node->Input(0);
  node->RedirectChildren(parent,node->InputConnNum(0));
  node->m_poss->DeleteChildAndCleanUp(node);
}


FindMidDistributions::FindMidDistributions(DistType srcType, DistType midType, DistType destType) 
  : m_srcType(srcType), m_midType(midType), m_destType(destType)
{
  m_type = "Finding middle distribution " + DistTypeToStr(m_srcType) + " -> " +
    DistTypeToStr(m_midType) + " -> " + DistTypeToStr(m_destType);
}


DLANode* FindRedistribution(DLANode *root, unsigned int num, DLANode *ignore, DistType type, bool goUp, unsigned int &outNum)
{
  if (!root)
    return NULL;
  if (root->GetDistType(num) == type) {
    outNum = num;
    return root;
  }
  NodeConnVecConstIter iter = root->m_children.begin();
  for( ; iter != root->m_children.end(); ++iter) {
    if ((*iter)->m_num == num) {
      DLANode *output = (DLANode*)(*iter)->m_n;
      if (output != ignore) {
        if ((output->GetNodeClass() == RedistNode::GetClass())) {
          if (output->GetDistType(0) == type) {
            outNum = 0;
            return output;
          }
          DLANode *find = FindRedistribution(output, 0, ignore, type, false, outNum);
          if (find) {
            return find;
          }
        }
      }
    }
  }
  if (goUp && root->GetNodeClass() == RedistNode::GetClass())
    return FindRedistribution((DLANode*)(root->Input(0)), root->InputConnNum(0), root, type, true, outNum);
  else
    return NULL;
}


bool FindMidDistributions::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    return false;
  DLANode *ddla = (DLANode*)node;
  if (ddla->GetDistType(0) != m_destType)
    return false;
  DLANode *parent = (DLANode*)(ddla->Input(0));
  if (!parent
      || ((parent)->GetDistType(node->InputConnNum(0)) != m_srcType))
    return false;
  for(unsigned int i = 0; i < ddla->m_children.size(); ++i) {
    if (ddla->Child(i)->GetNodeClass() == RedistNode::GetClass()) {
      RedistNode *tmp = (RedistNode*)(ddla->Child(i));
      if (tmp->m_destType == m_midType) {
	return true;
      }
    }
  }
  unsigned int num;
  return FindRedistribution(parent, ddla->InputConnNum(0), ddla, m_midType, true, num) != NULL;
}

void FindMidDistributions::Apply(Node *node) const
{
  DLANode *ddla = (DLANode*)node;
  DLANode *parent = (DLANode*)(ddla->Input(0));
  for(unsigned int i = 0; i < ddla->m_children.size(); ++i) {
    if (ddla->Child(i)->GetNodeClass() == RedistNode::GetClass()) {
      RedistNode *tmp = (RedistNode*)(ddla->Child(i));
      if (tmp->m_destType == m_midType) {
	tmp->ChangeInput2Way(ddla, ddla->ChildConnNum(i), parent, ddla->InputConnNum(0));
	ddla->ChangeInput2Way(parent, ddla->InputConnNum(0), tmp, 0);
	if (ddla->m_children.size() == 0)
	  ddla->m_poss->DeleteChildAndCleanUp(ddla);
	return;
      }
    }
  }
  unsigned int num;
  DLANode *newParent = FindRedistribution(parent, ddla->InputConnNum(0), ddla, m_midType, true, num);
  if (!newParent)
    throw;
  unsigned int oldNum = node->InputConnNum(0);
  node->ChangeInput1Way(parent, oldNum, newParent, num);
  parent->RemoveChild(node, oldNum);
  if (parent->m_children.empty()) {
    parent->m_poss->DeleteChildAndCleanUp(parent);
  }
}

ReplaceWithTrans::ReplaceWithTrans(DistType srcType, DistType origDestType, DistType newDestType)
  : m_srcType(srcType), m_origDestType(origDestType), m_newDestType(newDestType)
{
  m_type = "Replacing " + DistTypeToStr(m_srcType) + " -> " +
    DistTypeToStr(m_origDestType) + " with transpose to " + DistTypeToStr(m_newDestType);
}

bool ReplaceWithTrans::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    return false;
  DLANode *ddla = (DLANode*)node;
  if (ddla->GetDistType(0) != m_origDestType)
    return false;
  Node *parent = ddla->Input(0);
  if (!parent
      || (((DLANode*)parent)->GetDistType(node->InputConnNum(0)) != m_srcType))
    return false;
  
  return true;
}

void ReplaceWithTrans::Apply(Node *node) const
{
  Node *parent = node->Input(0);
  RedistNode *newNode = new RedistNode(m_newDestType);
  node->m_poss->AddNode(newNode);
  newNode->AddInput(parent,node->InputConnNum(0));
  node->RedirectChildren(newNode,0);
  node->m_poss->DeleteChildAndCleanUp(node);
}

void TransTransformation::Apply(Node *node) const
{
  PreApply(node);
  Node *in = node->Input(m_argNum);
  if (in->GetNodeClass() == RedistNode::GetClass()) {
    RedistNode *oldRedist = (RedistNode*)(node->Input(m_argNum));
    RedistNode *newRedist = oldRedist->CreateTrans(m_trans);
    newRedist->AddInput(oldRedist->Input(0),oldRedist->InputConnNum(0));
    node->ChangeInput1Way(oldRedist, 0, newRedist, 0);
    oldRedist->RemoveChild(node, node->InputConnNum(m_argNum));
    node->m_poss->AddNode(newRedist);
    if (oldRedist->m_children.empty())
      oldRedist->m_poss->DeleteChildAndCleanUp(oldRedist);
  }
  else if (in->GetNodeClass() == MakeTrapNode::GetClass()) {
    MakeTrapNode *trap = (MakeTrapNode*)(node->Input(m_argNum));
    if (trap->Input(0)->GetNodeClass() != RedistNode::GetClass())
      throw;
    RedistNode *oldRedist = (RedistNode*)(trap->Input(0));
    RedistNode *newRedist = oldRedist->CreateTrans(m_trans);
    newRedist->AddInput(oldRedist->Input(0),oldRedist->InputConnNum(0));
    trap->ChangeInput1Way(oldRedist, 0, newRedist, 0);
    oldRedist->RemoveChild(trap, trap->InputConnNum(0));
    node->m_poss->AddNode(newRedist);
    if (oldRedist->m_children.empty())
      oldRedist->m_poss->DeleteChildAndCleanUp(oldRedist);
  }
  else
    throw;
  PostApply(node);
  
}

string TransTransformation::GetType() const
{
  string str = GetTransType();
  str = str + " Trans " + (char)(m_argNum + 48);
  str = str + " " + TransToStr(m_trans);
  return str;
}


string RedistTrans::GetTransType() const
{
  return "redist";
}

bool RedistTrans::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    return false;
  RedistNode *redist = (RedistNode*)node;
  DLANode *source = (DLANode*)redist->Input(0);
  switch(redist->m_destType) {
  case(D_STAR_VR):
    return redist->InputDistType(0) == D_STAR_MR && source->CanTrans();
  case(D_VR_STAR):
    return redist->InputDistType(0) == D_MR_STAR && source->CanTrans();
  case(D_STAR_VC):
    return redist->InputDistType(0) == D_STAR_MC && source->CanTrans();
  case(D_VC_STAR):
    return redist->InputDistType(0) == D_MC_STAR && source->CanTrans();
  default:
    return false;
  }
}

bool RedistTrans::WorthApplying(const Node *node) const
{
  return true;
}

RedistNode::RedistNode(DistType destType)
  : m_destType(destType), m_mSizes(NULL), m_nSizes(NULL)
{
}

RedistNode::~RedistNode()
{
  if (m_mSizes) {
    delete m_mSizes;
    m_mSizes = NULL;
    delete m_nSizes;
    m_nSizes = NULL;
  }
}

void RedistNode::Duplicate(const Node *orig, bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig,shallow, possMerging);
  const RedistNode *origNode = (RedistNode*)orig;
  m_destType = origNode->m_destType;
}

NodeType RedistNode::GetType() const
{
  if (m_inputs.empty()) {
    return "RedistNode to " + DistTypeToStr(m_destType) +
      " without parent";
  }
  else {
    Node *parent = Input(0);
    DistType type = ((DLANode*)parent)->GetDistType(InputConnNum(0));
    return  "RedistNode " + DistTypeToStr(type) +
      " -> " + DistTypeToStr(m_destType);
  }
}

void RedistNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();
    if (m_inputs.size() != 1) {
      cout << "m_inputs.size() != 1\n";
      throw;
    }
  
    if (!m_children.size())
      throw;

    DLANode *parent = (DLANode*)Input(0);
    DistType m_srcType = parent->GetDistType(InputConnNum(0));
    const Sizes *m = GetM(0);
    const Sizes *n = GetN(0);
    parent->Prop();
    m_cost = GetCost(m_srcType, m_destType, m, n);
  }

}

Cost RedistNode::GetCost(DistType srcType, DistType destType, const Sizes *m, const Sizes *n)
{
  Cost cost = -1;
  try {
      switch (destType) {
      case(D_MR_MC_H):
	cost = GetRedistCost<D_MR_MC_H>(srcType, m, n);
	break;
      case (D_MR_MC_T):
	cost = GetRedistCost<D_MR_MC_T>(srcType, m, n);
	break;
      case(D_MC_STAR):
	cost = GetRedistCost<D_MC_STAR>(srcType, m, n);
	break;
      case(D_VC_STAR):
	cost = GetRedistCost<D_VC_STAR>(srcType, m, n);
	break;
      case(D_STAR_VC):
	cost = GetRedistCost<D_STAR_VC>(srcType, m, n);
	break;
      case (D_MC_MR):
	cost = GetRedistCost<D_MC_MR>(srcType, m, n);
	break;
      case (D_VR_STAR):
	cost = GetRedistCost<D_VR_STAR>(srcType, m, n);
	break;
      case (D_STAR_VR):
	cost = GetRedistCost<D_STAR_VR>(srcType, m, n);
	break;
      case (D_MR_STAR):
	cost = GetRedistCost<D_MR_STAR>(srcType, m, n);
	break;
      case (D_MR_STAR_T):
	cost = GetRedistCost<D_MR_STAR_T>(srcType, m, n);
	break;
      case (D_MR_STAR_H):
	cost = GetRedistCost<D_MR_STAR_H>(srcType, m, n);
	break;
      case (D_STAR_STAR):
	cost = GetRedistCost<D_STAR_STAR>(srcType, m, n);
	break;
      case (D_STAR_MC):
	cost = GetRedistCost<D_STAR_MC>(srcType, m, n);
	break;
      case (D_STAR_MC_T):
	cost = GetRedistCost<D_STAR_MC_T>(srcType, m, n);
	break;
      case (D_STAR_MC_H):
	cost = GetRedistCost<D_STAR_MC_H>(srcType, m, n);
	break;
      case (D_STAR_MR):
	cost = GetRedistCost<D_STAR_MR>(srcType, m, n);
	break;
      case (D_STAR_MR_T):
	cost = GetRedistCost<D_STAR_MR_T>(srcType, m, n);
	break;
      case (D_STAR_MR_H):
	cost = GetRedistCost<D_STAR_MR_H>(srcType, m, n);
	break;
      case (D_MR_MC):
	cost = GetRedistCost<D_MR_MC>(srcType, m, n);
	break;
      case (D_VC_STAR_H):
	cost = GetRedistCost<D_VC_STAR_H>(srcType, m, n);
	break;
      case (D_VC_STAR_T):
	cost = GetRedistCost<D_VC_STAR_T>(srcType, m, n);
	break;
      default:
	cout << "Not handling " << DistTypeToStr(srcType) << " -> " << DistTypeToStr(destType) << endl;
	throw;
	break;
      }
    }
    catch (...) {
      throw;
    }
  
  if(m && n && !cost) {
    cout << DistTypeToStr(srcType) << " to " << DistTypeToStr(destType) << " cost = " << cost << endl;
    cout << m << " by " << n << endl;
    cout << "blah!\n";
  }
  
  if(!IsValidCost(cost)) {
    cout << cost << endl;
    throw;
  }
  return cost;
}

const Sizes* RedistNode::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputM(0);
}

const Sizes* RedistNode::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputN(0);
}

Name RedistNode::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  Name name = GetInputName(0);
  name.m_type = m_destType;
  return name;
}

void RedistNode::PrintCode(IndStream &out)
{  
  out.Indent();
  DistType in = InputDistType(0);
  if ((in == D_MR_MC_H && m_destType == D_MC_MR) || 
      (in == D_VC_STAR_H && m_destType == D_STAR_MC)) {
    *out << "Adjoint( " << GetInputNameStr(0) << ".LockedLocalMatrix(), " 
	 << GetName(0).str() << " );\n";
  }
  else if ((in == D_MR_MC_T && m_destType == D_MC_MR) ||
	   (in == D_VC_STAR_T && m_destType == D_STAR_MC)) {
    *out << "Transpose( " << GetInputNameStr(0) << ".LockedLocalMatrix(), " 
    << GetName(0).str() << ".LocalMatrix() );\n";
  }
  else {
    switch(m_destType) {
    case (D_MC_STAR_T) :
    case (D_STAR_MC_T) :
    case (D_MR_STAR_T) :
    case (D_STAR_MR_T) :
      *out << GetName(0).str()
	   << ".TransposeFrom( "
	   << GetInputName(0).str()
	   << " );\n";
      break;
    case (D_MC_STAR_H) :
    case (D_STAR_MC_H) :
    case (D_MR_STAR_H) :
    case (D_STAR_MR_H) :
      *out << GetName(0).str()
	   << ".AdjointFrom( "
	   << GetInputName(0).str()
	   << " );\n";
      break;
    default:
      switch(((DLANode*)Input(0))->GetDistType(InputConnNum(0))) {
      case (D_MC_STAR_T) :
      case (D_STAR_MC_T) :
      case (D_MR_STAR_T) :
      case (D_STAR_MR_T) :
	*out << GetName(0).str()
	     << ".TransposeFrom( "
	     << GetInputName(0).str()
	     << " );\n";
	break;
      case (D_MC_STAR_H) :
      case (D_STAR_MC_H) :
      case (D_MR_STAR_H) :
      case (D_STAR_MR_H) :
	*out << GetName(0).str()
	     << ".AdjointFrom( "
	     << GetInputName(0).str()
	     << " );\n";
	break;
      default:
	*out << GetName(0).str()
	     << " = "
	     << GetInputName(0).str()
	     << ";\n";
	break;
      }
    }
  }
}

bool RedistNode::CanTrans() const
{
  return ::CanTrans(InputDistType(0),m_destType);
}

RedistNode* RedistNode::CreateTrans(Trans trans)
{
  if (trans != TRANS && trans != CONJTRANS) {
    cout << "Bad transpose!\n";
    return NULL;
  }
  if (!CanTrans()) {
    cout << "Tried to transpose when it wasn't possible\n";
    cout << DistTypeToStr(InputDistType(0)) << endl;
    cout << DistTypeToStr(m_destType) << endl;
    throw;
  }
  switch (m_destType) {
    case (D_STAR_MR):
      return new RedistNode(trans==TRANS ? D_MR_STAR_T : D_MR_STAR_H);
    case (D_STAR_MC):
      return new RedistNode(trans==TRANS ? D_MC_STAR_T : D_MC_STAR_H);
    case (D_MC_STAR):
      return new RedistNode(trans==TRANS ? D_STAR_MC_T : D_STAR_MC_H);
    case (D_MR_STAR):
      return new RedistNode(trans==TRANS ? D_STAR_MR_T : D_STAR_MR_H);
    default:
      cout << "Huh?";
      return NULL;
  }
}


void RedistNode::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_destType);
}


void RedistNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  READ(m_destType);
}

void SumScatterNode::Duplicate(const Node *orig,bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
  const SumScatterNode *origNode = (SumScatterNode*)orig;
  m_coeff = origNode->m_coeff;
}

bool SumScatterNode::KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const
{
  if (Input(1) == input && InputConnNum(1) == numIn) {
    numOut = 0;
    return true;
  }
  else
    return false;
}

bool SumScatterNode::Overwrites(const Node *input, unsigned int num) const
{
  const NodeConn *conn = m_inputs[1];
  return conn->m_n == input && conn->m_num == num;
}

NodeType SumScatterNode::GetType() const
{
  return "SumScatterNode";
}

void SumScatterNode::Prop()
{
  if (!IsValidCost(m_cost)) {
    DLANode::Prop();
    
  if (m_inputs.size() != 2)
    cout << "m_inputs.size() != 2\n";
    

    DLANode *from = (DLANode*)Input(0);
    DLANode *to = (DLANode*)Input(1);
  

    from->Prop();
    to->Prop();
    
    DistType destType = to->GetDistType(InputConnNum(1));
    
    DistType m_srcType = from->GetDistType(InputConnNum(0));

    const Sizes *localMs = LocalM(0);
    const Sizes *localNs = LocalN(0);
    
    m_cost = GetCost(localMs, localNs, destType, m_srcType);
  }
}

Cost SumScatterNode::GetCost(const Sizes *localMs, const Sizes *localNs, DistType destType, DistType srcType)
{
  Cost cost = 0;
  unsigned int length = localMs->NumSizes();

  if (length != localNs->NumSizes())
    throw;

  for(unsigned int i = 0; i < length; ++i) {
    Size localM = (*localMs)[i];
    Size localN = (*localNs)[i];
    
    if (srcType == D_MC_STAR && destType == D_MC_MR) {
      cost += ( C * CopyCost(localM, localN, ONE, ONE) )
	+ ReduceScatter(C * localM * localN, C)
	+ CopyCost( localM, localN, ONE, ONE)
	+ (TWO * GAMMA * localM * localN);
    }
    else if (srcType == D_MR_STAR && destType == D_MR_MC) {
      cost += ( R * CopyCost(localM, localN, ONE, ONE) )
	+ ReduceScatter(R * localM * localN, R)
	+ CopyCost( localM, localN, ONE, ONE)
	+ (TWO * GAMMA * localM * localN);
    }
    else if (srcType == D_STAR_STAR && 
	     (destType == D_MC_MR || destType == D_MR_MC)) 
      {
	cost += ( C * R * CopyCost( localM, localN, R, ONE))
	  + ReduceScatter(P * localM * localN, P)
	  + CopyCost( localM, localN, ONE, ONE )
	  + (TWO * GAMMA * localM * localN);
      }
    else if (srcType == D_STAR_MR && destType == D_MC_MR) {
      cost += (R * CopyCost(localM, localN, R, ONE))
	+ ReduceScatter(R*localM*localN, R) 
	+ CopyCost( localM, localN, ONE, ONE )
	+ (TWO * GAMMA * localM * localN);
    }
    else if (srcType == D_MR_STAR_T && destType == D_MC_MR) {
      cost += ( R * CopyCost(localN, localM, ONE, ONE) )
	+ ReduceScatter(R * localM * localN, R)
	+ CopyCost( localN, localM, ONE, ONE)
	+ CopyCost( localN, localM, ONE, localM );
    }
    else if (srcType == D_MR_STAR_H && destType == D_MC_MR) {
      cost += ( R * CopyCost(localN, localM, ONE, ONE) )
	+ ReduceScatter(R * localM * localN, R)
	+ CopyCost( localN, localM, ONE, ONE)
	+ CopyCost( localN, localM, ONE, localM );
    }   
    else if (srcType == D_STAR_MC && destType == D_MR_MC) {
      cost += (C * CopyCost(localM, localN, C, ONE))
	+ ReduceScatter(C*localM*localN, C) 
	+ CopyCost( localM, localN, ONE, ONE )
	+ (TWO * GAMMA * localM * localN);
    }
    else if ((srcType == D_MR_STAR_H && destType == D_MR_MC_H) ||
	     (srcType == D_MR_STAR_T && destType == D_MR_MC_T)) {
      cost += ( R * CopyCost(localM, localN, ONE, ONE) )
	+ ReduceScatter(R * localM * localN, R)
	+ CopyCost( localM, localN, ONE, ONE)
	+ (TWO * GAMMA * localM * localN);
    }
    else
      throw;
  }
  return cost;
}

const Sizes* SumScatterNode::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputM(1);
}

const Sizes* SumScatterNode::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputN(1);
}

const Sizes* SumScatterNode::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputLocalM(1);
}

const Sizes* SumScatterNode::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputLocalN(1);
}

Name SumScatterNode::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  DistType m_destType = ((DLANode*)Input(1))->GetDistType(InputConnNum(1));
  Name tmp = GetInputName(1);
  tmp.m_type = m_destType;
  return tmp;
}

void SumScatterNode::PrintCode(IndStream &out)
{  
  DistType inType = InputDistType(0);
  out.Indent();
  Trans trans = DistTransType(inType);
  if (trans != NORMAL) {
    Trans trans2 = DistTransType(InputDistType(1));
    if (trans2 != NORMAL) {
      if (trans != trans2)
	throw;
      trans = NORMAL;
    }
  }

  switch(trans) {
  case (NORMAL):
    *out << GetName(0).str()
	 << ".SumScatterUpdate( ";
    out << m_coeff;
    *out << ", " << GetInputName(0).str() << " );" << endl;
    break;
  case (TRANS):
    *out << GetName(0).str()
	 << ".TransposeSumScatterUpdate( ";
    out<< m_coeff;
    *out << ", " << GetInputName(0).str() << " );" << endl;
    break;
  case (CONJTRANS):
    *out << GetName(0).str()
	 << ".AdjointSumScatterUpdate( ";
    out << m_coeff;
    *out << ", " << GetInputName(0).str() << " );" << endl;
    break;
  default:
    throw;
  }
}

void SumScatterNode::FlattenCore(ofstream &out) const
{
  DLANode::FlattenCore(out);
  WRITE(m_coeff);
}

void SumScatterNode::UnflattenCore(ifstream &in, SaveInfo &info)
{
  DLANode::UnflattenCore(in, info);
  READ(m_coeff);
}

void SumScatterFrom::Prop()
{
  DLANode *from = (DLANode*)Input(0);
  DLANode *to = (DLANode*)Input(1);
  
  if (!IsValidCost(m_cost)) {
    if (m_inputs.size() != 2)
      cout << "m_inputs.size() != 2\n";

    DLANode::Prop();

    from->Prop();
    to->Prop();
    
    DistType destType = to->GetDistType(InputConnNum(1));
    DistType m_srcType = from->GetDistType(InputConnNum(0));

    const Sizes *localMs = LocalM(0);
    const Sizes *localNs = LocalN(0);

    m_cost = GetCost(destType, m_srcType, localMs, localNs);
  }
}

bool SumScatterFrom::Overwrites(const Node *input, unsigned int num) const
{
  const NodeConn *conn = m_inputs[1];
  return conn->m_n == input && conn->m_num == num;
}

bool SumScatterFrom::KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const
{
  if (Input(1) == input && InputConnNum(1) == numIn) {
    numOut = 0;
    return true;
  }
  else
    return false;
}

Cost SumScatterFrom::GetCost(DistType destType, DistType srcType, const Sizes *localMs, const Sizes *localNs)
{
  unsigned int length = localMs->NumSizes();
  if (length != localNs->NumSizes())
    throw;

  Cost cost = 0;

  for(unsigned int i = 0; i < length; ++i) {
    Size localM = (*localMs)[i];      
    Size localN = (*localNs)[i];
    
    if (srcType == D_STAR_MC && destType == D_MR_MC) {
      cost += (C * CopyCost(localM,localN,C, ONE))
	+ ReduceScatter(C * localM * localN, C)
	+ CopyCost( localM, localN, ONE, ONE);
    }
    else if (srcType == D_MR_STAR && destType == D_MR_MC) {
      cost += ( R * CopyCost(localM, localN, ONE, ONE) )
	+ ReduceScatter(R * localM * localN, R)
	+ CopyCost( localM, localN, ONE, ONE);
    }
    else if ((srcType == D_MC_STAR_H && destType == D_MC_MR_H) ||
	     (srcType == D_MC_STAR_T && destType == D_MC_MR_T)) {
      cost += (C * CopyCost(localM,localN,C, ONE))
	+ ReduceScatter(C * localM * localN, C)
	+ CopyCost( localM, localN, ONE, ONE);
    }
    else if (srcType == D_MC_STAR && destType == D_MC_MR) {
      cost += (C * CopyCost(localM,localN,ONE,ONE))
	+ ReduceScatter(C * localM * localN, C)
	+ CopyCost(localM, localN, ONE, ONE);
    }
    else if ((srcType == D_MR_STAR_T && destType == D_MR_MC_T) ||
	     (srcType == D_MR_STAR_H && destType == D_MR_MC_H)) {
      cost += ( R * CopyCost(localM, localN, ONE, ONE) )
	+ ReduceScatter(R * localM * localN, R)
	+ CopyCost( localM, localN, ONE, ONE);
    }
    else
      throw;
  }

  return cost;
}

const Sizes* SumScatterFrom::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputM(1);
}

const Sizes* SumScatterFrom::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputN(1);
}

const Sizes* SumScatterFrom::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputLocalM(1);
}

const Sizes* SumScatterFrom::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputLocalN(1);
}

Name SumScatterFrom::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  DistType m_destType = ((DLANode*)Input(1))->GetDistType(InputConnNum(1));
  Name tmp = GetInputName(1);
  tmp.m_type = m_destType;
  return tmp;
}

void SumScatterFrom::PrintCode(IndStream &out)
{  
  out.Indent();
  *out << GetName(0).str()
  << ".SumScatterFrom( "
  << GetInputName(0).str() << " );" << endl;
}

void SumOverCommNode::Duplicate(const Node *orig,bool shallow, bool possMerging)
{
  DLANode::Duplicate(orig, shallow, possMerging);
}

bool SumOverCommNode::KeepsInputVarLive(Node *input, unsigned int numIn, unsigned int &numOut) const
{
  if (Input(0) == input && InputConnNum(0) == numIn) {
    numOut = numIn;
    return true;
  }
  else
    throw;
}

bool SumOverCommNode::Overwrites(const Node *input, unsigned int num) const
{
  const NodeConn *conn = m_inputs[0];
  return conn->m_n == input && conn->m_num == num;
}

NodeType SumOverCommNode::GetType() const
{
  return "SumOverCommNode";
}

void SumOverCommNode::Prop()
{
  if (!IsValidCost(m_cost)) {
     DLANode::Prop();
    if (m_inputs.size() != 1) {
      cout << "m_inputs.size() != 1\n";
      throw;
    }   
 
   DistType destType = InputDistType(0);
    if (GetRowDist(destType) != STAR && GetColDist(destType) != STAR)
      throw;

    const Sizes *localMs = LocalM(0);
    const Sizes *localNs = LocalN(0);
    
    m_cost = GetCost(destType, localMs, localNs);
  }
}

Cost SumOverCommNode::GetCost(DistType destType, const Sizes *localMs, const Sizes *localNs)
{
  unsigned int length = localMs->NumSizes();
  if (length != localNs->NumSizes())
    throw;

  Cost cost = 0;

  for(unsigned int i = 0; i < length; ++i) {
    Size localM = (*localMs)[i];
    Size localN = (*localNs)[i];
    
    if (destType == D_MC_STAR) {
      cost += CopyCost(localM, localN, ONE, ONE) 
	+ AllReduce(localM * localN, C)
	+ CopyCost( localM, localN, ONE, ONE);
    }
    else if (destType == D_MR_STAR) {
      cost += CopyCost(localM, localN, ONE, ONE) 
	+ AllReduce(localM * localN, R)
	+ CopyCost( localM, localN, ONE, ONE);
    }
    else
      throw;
  }
  return cost;
}

const Sizes* SumOverCommNode::GetM(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputM(0);
}

const Sizes* SumOverCommNode::GetN(unsigned int num) const
{
  if (num > 0)
    throw;
  return GetInputN(0);
}



const Sizes* SumOverCommNode::LocalM(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputLocalM(0);
}

const Sizes* SumOverCommNode::LocalN(unsigned int num) const
{
  if (num > 0)
    throw;
  return InputLocalN(0);
}

Name SumOverCommNode::GetName(unsigned int num) const
{
  if (num > 0)
    throw;
  DistType m_destType = ((DLANode*)Input(0))->GetDistType(InputConnNum(0));
  Name tmp = GetInputName(0);
  tmp.m_type = m_destType;
  return tmp;
}

void SumOverCommNode::PrintCode(IndStream &out)
{
  DistType type = InputDistType(0);
  out.Indent();
  *out << GetName(0).str();
  if (type == D_MC_STAR || type == D_STAR_MC)
    *out << ".SumOverRow();\n";
  else if (type == D_MR_STAR || type == D_STAR_MR)
    *out << ".SumOverCol();\n";
  else 
    throw;
}

void RedistNode::ClearDataTypeCache()
{
  if (m_mSizes) {
    delete m_mSizes;
    m_mSizes = NULL;
    delete m_nSizes;
    m_nSizes = NULL;
  }
}

void RedistNode::BuildDataTypeCache()
{
  if (m_mSizes)
    return;
  DLANode *in = (DLANode*)Input(0);
  unsigned int num = InputConnNum(0);
  const Sizes *ms = in->GetM(num);
  const Sizes *ns = in->GetN(num);
  if (ms->NumSizes() != ns->NumSizes()) {
    cout << ms->NumSizes() << endl;
    cout << ns->NumSizes() << endl;
    throw;
  }
  m_mSizes = new Sizes;
  m_nSizes = new Sizes;
  GetLocalSizes(m_destType, ms, ns, *m_mSizes, *m_nSizes);
}

const Sizes* RedistNode::LocalM(unsigned int num) const
{
  if (m_mSizes)
    return m_mSizes;
  else {
    cout << "no size cache on RedistNode\n";
    throw;
  }
}

const Sizes* RedistNode::LocalN(unsigned int num) const
{
  if (m_nSizes)
    return m_nSizes;
  else {
    cout << "no size cache on RedistNode\n";
    throw;
  }
}




bool CanUseUp(const Node *root)
{
  if (root->GetNodeClass() != RedistNode::GetClass()) {
    cout << "Root isn't redist node\n";
    throw;
    //    return false;
  }
  else {
    if (root->Input(0)->GetNodeClass() == RedistNode::GetClass())
      return true;
    const DLANode *prev = (DLANode*)root;
    const DLANode *curr = (DLANode*)root;
    unsigned int outputNum = root->InputConnNum(0);
    while (curr->GetNodeClass() == RedistNode::GetClass()) {
      for (unsigned int i = 0; i < curr->m_children.size(); ++i) {
        const Node *child = curr->Child(i);
        if (child != prev && child->GetNodeClass() == RedistNode::GetClass())
          return true;
      }
      prev = curr;
      if(curr->m_inputs.size() <= 0)
        return false;
      outputNum = curr->InputConnNum(0);
      curr = (DLANode*)curr->Input(0);
      if (!curr)
        throw;
    }
    return curr->NumChildrenOfOutput(outputNum) > 1;
  }
}

bool UniqueTransTrans::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass()) 
    return false;
  RedistNode *redist = (RedistNode*)node;
  if (redist->m_destType == D_STAR_MC) {
    const DLANode *par = (DLANode*)(redist->Input(0));
    unsigned int num = redist->InputConnNum(0);
    NodeConnVecConstIter iter = par->m_children.begin();
    for(; iter != par->m_children.end(); ++iter) {
      if ((*iter)->m_num == num) {
	const DLANode *child = (DLANode*)(*iter)->m_n;
	if (child->GetNodeClass() == RedistNode::GetClass()) {
	  const RedistNode *childRedist = (RedistNode*)child;
	  if (childRedist->m_destType == D_MR_STAR_H || childRedist->m_destType == D_MR_STAR_T)
	    return true;
	}
      }
    }
  }
  return false;
}

void UniqueTransTrans::Apply(Node *node) const 
{
  if (node->GetNodeClass() != RedistNode::GetClass()) 
    throw;
  RedistNode *redist = (RedistNode*)node;
  if (redist->m_destType == D_STAR_MC) {
    DLANode *par = (DLANode*)(redist->Input(0));
    unsigned int num = redist->InputConnNum(0);
    NodeConnVecIter iter = par->m_children.begin();
    for(; iter != par->m_children.end(); ++iter) {
      if ((*iter)->m_num == num) {
	DLANode *child = (DLANode*)(*iter)->m_n;
	if (child->GetNodeClass() == RedistNode::GetClass()) {
	  RedistNode *childRedist = (RedistNode*)child;
	  if (childRedist->m_destType == D_MR_STAR_H || childRedist->m_destType == D_MR_STAR_T) {
	    RedistNode *redist 
	      = new RedistNode(childRedist->m_destType == D_MR_STAR_H ? D_VC_STAR_H : D_VC_STAR_T);
	    redist->AddInput(childRedist, 0);
	    node->ChangeInput2Way(par, num, redist, 0);
	    node->m_poss->AddNode(redist);
	    return;
	  }
	}
      }
    }
  }
  throw;
}




UseTransposedRedist::UseTransposedRedist(DistType destType, DistType srcType1, DistType srcType2)
{
  m_destType = destType;
  m_srcType1 = srcType1;
  m_srcType2 = srcType2;
}

string UseTransposedRedist::GetType() const 
{ 
  return "Use transposed redist for " + DistTypeToStr(m_destType); 
}

bool UseTransposedRedist::CanApply(const Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    return false;
  RedistNode *redist = (RedistNode*)node;
  if (redist->m_destType != m_destType)
    return false;
  if (redist->InputDistType(0) == m_srcType1 ||
      redist->InputDistType(0) == m_srcType2)
    return false;
  DLANode *parent = (DLANode*)(redist->Input(0));
  unsigned int outNum;
  if (FindRedistribution(parent, redist->InputConnNum(0), redist,
			 m_srcType1, true, outNum)) {
    return true;
  }
  if (FindRedistribution(parent, redist->InputConnNum(0), redist,
			 m_srcType2, true, outNum)) {
    return true;
  }
  return false;
}

void UseTransposedRedist::Apply(Node *node) const
{
  if (node->GetNodeClass() != RedistNode::GetClass())
    throw;
  RedistNode *redist = (RedistNode*)node;
  if (redist->m_destType != m_destType)
    throw;
  DLANode *parent = (DLANode*)(redist->Input(0));
  unsigned int outNum;
  DLANode *newParent = NULL;
  newParent = FindRedistribution(parent, redist->InputConnNum(0), redist,
			    m_srcType1, true, outNum);
  if (!newParent)
    newParent = FindRedistribution(parent, redist->InputConnNum(0), redist,
			      m_srcType2, true, outNum);
  if (!newParent)
    throw;
  redist->ChangeInput2Way(parent, redist->InputConnNum(0), newParent, outNum);
  if (parent->m_children.empty()) {
    parent->m_poss->DeleteChildAndCleanUp(parent);
  }
}


template class ExpandRedistribution<D_MC_MR, D_VR_STAR>;
template class ExpandRedistribution<D_MC_MR, D_MC_STAR>;
template class ExpandRedistribution<D_MC_MR, D_STAR_MR>;
template class ExpandRedistribution<D_MC_MR, D_MR_STAR>;
template class ExpandRedistribution<D_MC_MR, D_STAR_MC>;
template class ExpandRedistribution<D_VC_STAR,D_MR_STAR>;

#endif
