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

#include "sizes.h"
#include "costs.h"

SizesCache::~SizesCache()
{
  for(auto& elem : m_midSizesMap) {
    elem.second->m_cached = false;
    delete elem.second;
  }
  for(auto& elem : m_rangeMap) {
    elem.second->m_cached = false;
    delete elem.second;
  }

  for(auto& elem : m_repSizesMap) {
    elem.second->m_cached = false;
    delete elem.second;
  }

  for(auto& elem : m_numItersMap) {
    delete elem.second;
  }

  for(auto& elem : m_constMap) {
    elem.second->m_cached = false;
    delete elem.second;
  }

  for(auto& elem : m_otherSizes) {
    delete elem;
  }
}

void SizesCache::TakeSize(SizeList *size)
{
  size->SetCached();
  m_otherSizes.push_back(size);
}

const SizeList* SizesCache::GetCachedMidSize(const SizeList *parent,
					     Size size)
{
  SizesT<Size> val;
  val.parent = parent;
  val.size = size;
  SizesSizeMapIter find = m_midSizesMap.find(val);
  if (find != m_midSizesMap.end())
    return find->second;
  else {
    if (!parent->IsCached())
      throw;
    SizeList *sizes = new SizeList;
    unsigned int numExecs = parent->NumSizes();

    for(unsigned int i = 0; i < numExecs; ++i) {
      Size len = (*parent)[i];
      sizes->AddMidSizes(size, len);
    }
    sizes->SetCached();
    m_midSizesMap[val] = sizes;
    return sizes;
  }
}

const SizeList* SizesCache::GetCachedRange(bool start,
					   const SizeList *parent,
					   int stride)
{
  if (!start)
    stride *= -1;
  SizesT<int> val;
  val.parent = parent;
  val.size = stride;
  SizesIntMapIter find = m_rangeMap.find(val);
  if (find != m_rangeMap.end())
    return find->second;
  else {
    if (!parent->IsCached())
      throw;

    SizeList *sizes = new SizeList;
    unsigned int numExecs = parent->NumSizes();

    for(unsigned int i = 0; i < numExecs; ++i) {
      Size len = (*parent)[i];
      if (start)
	sizes->AddSizesWithLimit(0,stride,len);
      else
	sizes->AddSizesWithLimit(len,stride,0);
    }

    sizes->SetCached();

    m_rangeMap[val] = sizes;

    return sizes;
  }
}

const SizeList* SizesCache::GetCachedRepeatedSize(const SizeList *parent,
						  const SizeList *controlParent,
						  const int size)
{
  RepSizesData val;
  val.parent = parent;
  val.controlParent = controlParent;
  val.size = size;

  RepSizesMapIter find = m_repSizesMap.find(val);
  if (find != m_repSizesMap.end())
    return find->second;
  else {
    if (!parent->IsCached())
      throw;
    if (!controlParent->IsCached())
      throw;
    SizeList *sizes = new SizeList;

    const NumItersVec *numIters = GetNumItersVec(controlParent, size);
    
    unsigned int numExecs = parent->NumSizes();
    if (numExecs != numIters->size())
      throw;

    for(unsigned int i = 0; i < numExecs; ++i) {
      Size len = (*parent)[i];
      sizes->AddRepeatedSizes(len, (*numIters)[i]);
    }
    
    sizes->SetCached();
    m_repSizesMap[val] = sizes;
    return sizes;
  }
}

const SizeList* SizesCache::GetConstSize(Size size)
{
  SizeMapIter find = m_constMap.find(size);
  if (find != m_constMap.end())
    return find->second;
  else {
    SizeList *sizes = new SizeList;
    sizes->AddRepeatedSizes(size, 1);
    sizes->SetCached();
    m_constMap[size] = sizes;
    return sizes;
  }
}


#if DOTENSORS
const SizeList* SizesCache::GetCachedDistSize(const SizeList *parent,
					      DistEntry entry)
{
  DistSizesData val;
  val.parent = parent;
  val.entry = entry;

  DistSizesMapIter find = m_distSizesMap.find(val);
  if (find != m_distSizesMap.end())
    return find->second;
  else {
    SizeList *size = new SizeList;
    *size = *parent;
    if (!entry.IsStar()) {
      DimVec vec = entry.DistEntryDims();
      unsigned int coef = 1;
      DimVecIter iter = vec.begin();
      for(; iter != vec.end(); ++iter) {
	coef *= GridLens[*iter];
      }
      size->SetCoeff(1.0 / coef);
    }
    m_distSizesMap[val] = size;
    size->SetCached();
    return size;
  }
}
#endif


const NumItersVec* SizesCache::GetNumItersVec(const SizeList *controlParent,
					      const int size)
{
  SizesT<int> val;
  val.parent = controlParent;
  val.size = size;
  SizesNumIterMapIter iter = m_numItersMap.find(val);
  if (iter != m_numItersMap.end())
    return iter->second;
  else {
    NumItersVec *vec = new NumItersVec;
    unsigned int numExecs = controlParent->NumSizes();
    for(unsigned int i=0; i < numExecs; ++i) {
      vec->push_back((unsigned int)ceil((*controlParent)[i] / (double)size));
    }
    m_numItersMap[val] = vec;
    return vec;
  }
}


const SizeList* SizesCache::GetCachedRepeatedSize(Size size,
						  unsigned int numRepeats)
{
  RepeatedData val;
  val.totSize = size;
  val.reps = numRepeats;

  RepeatedMapIter find = m_constRepMap.find(val);
  if (find != m_constRepMap.end())
    return find->second;
  else {
    SizeList *newSize = new SizeList;
    newSize->AddRepeatedSizes(size, numRepeats);
    newSize->SetCached();
    m_constRepMap[val] = newSize;
    return newSize;
  }

}


const SizeList* SizesCache::GetCachedSplitSize(bool start,
					       const SizeList *parent,
					       int splitFactor)
{
  SizesT<int> val;
  val.parent = parent;
  val.size = splitFactor;
  if (start) {
    SizesIntMapIter find = m_splitMapStart.find(val);
    if (find != m_splitMapStart.end())
      return find->second;
  }
  else {
    SizesIntMapIter find = m_splitMapEnd.find(val);
    if (find != m_splitMapEnd.end())
      return find->second;
  }

  SizeList *startSizes = new SizeList();
  SizeList *endSizes = new SizeList();

  for (auto entry : parent->m_entries) {
    if (entry->m_type != REPEATEDSIZES)
      LOG_FAIL("havne't implement other splitting code");
    int numIterations = entry->NumSizesPerRepeat();
    Size sizeOfEachIteration = entry->m_valA;
    Size sizeOfStartIterations = floor(sizeOfEachIteration / (double)splitFactor) * splitFactor;
    Size sizeOfEndIterations = sizeOfEachIteration - sizeOfStartIterations;

    auto startEnt = new SizeEntry();
    startEnt->SetRepeatedSizes(sizeOfStartIterations, numIterations);
    startEnt->m_repeats = entry->m_repeats;
    startSizes->m_entries.push_back(startEnt);

    auto endEnt = new SizeEntry();
    endEnt->SetRepeatedSizes(sizeOfEndIterations, numIterations);
    endEnt->m_repeats = entry->m_repeats;
    endSizes->m_entries.push_back(endEnt);
  }

  m_splitMapStart[val] = startSizes;
  m_splitMapEnd[val] = endSizes;
  startSizes->SetCached();
  endSizes->SetCached();

  if (start)
    return startSizes;
  else
    return endSizes;
}
