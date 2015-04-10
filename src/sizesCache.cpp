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

#include "sizesCache.h"

SizesCache::~SizesCache()
{
  for(auto elem : m_midSizesMap)
    delete elem.second;
  for(auto elem : m_rangeStartMap)
    delete elem.second;
  for(auto elem : m_rangeEndMap)
    delete elem.second;
  for(auto elem : m_repSizesMap)
    delete elem.second;
  for(auto elem : m_numItersMap)
    delete elem.second;
}

const Sizes* SizesCache::GetCachedMidSize(const Sizes *parent,
					  Size size,
					  const NumItersVec *numIters)
{
  SizesT<Size> val;
  val.parent = parent;
  val.size = size;
  SizesSizeMapIter find = m_midSizesMap.find(val);
  if (find != m_midSizesMap.end())
    return find->second;
  else {
    Sizes *sizes = new Sizes;
    unsigned int numExecs = parent->NumSizes();
    if (numIters && numExecs != numIters->size())
      throw;
    for(unsigned int i = 0; i < numExecs; ++i) {
      Size len = (*parent)[i];
      unsigned int iters = sizes->AddMidSizes(size, len);
      if (numIters && (*numIters)[i] != iters) {
	throw;
      }
    }
    m_midSizesMap[val] = sizes;
    return sizes;
  }
}

const Sizes* SizesCache::GetCachedRange(bool start,
					const Sizes *parent,
					int stride,
					const NumItersVec *numIters)
{
  SizesT<int> val;
  val.parent = parent;
  val.size = stride;
  SizesIntMapIter find = (start ? 
			  m_rangeStartMap.find(val) :
			  m_rangeEndMap.find(val));
  if (start && find != m_rangeStartMap.end())
    return find->second;
  else if (!start && find != m_rangeEndMap.end())
    return find->second;
  else {
    Sizes *sizes = new Sizes;
    unsigned int numExecs = parent->NumSizes();

    NumItersVec *vec = NULL;
    if (m_numItersMap.find(val) == m_numItersMap.end()) {
      vec = new NumItersVec;
      vec->reserve(numExecs);
    }

    if (numIters && numExecs != numIters->size())
      throw;
    for(unsigned int i = 0; i < numExecs; ++i) {
      Size len = (*parent)[i];
      unsigned int iters;
      if (start) 
	iters = sizes->AddSizesWithLimit(len,-1*stride,0);
      else
	iters  = sizes->AddSizesWithLimit(0,stride,len);
      if (vec)
	vec->push_back(iters);
      if (numIters && (*numIters)[i] != iters) {
	throw;
      }
    }

    if (vec)
      m_numItersMap[val] = vec;

    if (start)
      m_rangeStartMap[val] = sizes;
    else
      m_rangeEndMap[val] = sizes;

    return sizes;
  }
}

const Sizes* SizesCache::GetCachedRepeatedSize(const Sizes *parent,
					       const Sizes *controlParent,
					       const Size size,
					       const NumItersVec *numIters)
{
  RepSizesData val;
  val.parent = parent;
  val.controlParent = controlParent;
  val.size = size;

  RepSizesMapIter find = m_repSizesMap.find(val);
  if (find != m_repSizesMap.end())
    return find->second;
  else {
    Sizes *sizes = new Sizes;
    
    unsigned int numExecs = parent->NumSizes();
    if (numExecs != numIters->size())
      throw;

    for(unsigned int i = 0; i < numExecs; ++i) {
      Size len = (*parent)[i];
      sizes->AddRepeatedSizes(len, (*numIters)[i]);
    }
    
    m_repSizesMap[val] = sizes;
    return sizes;
  }
}
