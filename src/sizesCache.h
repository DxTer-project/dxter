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

#include "sizes.h"
#include <map>
#include <vector>
#include "base.h"

using namespace std;

template <typename T>
struct SizesT {
  const SizeList *parent;
  T size;
};

template <typename T>
struct SizesTCompare {
  bool operator() (const SizesT<T>& lhs, const SizesT<T>& rhs) const{
    if (lhs.parent < rhs.parent) {
      return true;
    }
    else if (lhs.parent == rhs.parent) {
      return lhs.size < rhs.size;
    }
    else
      return false;
  }
};


typedef map<SizesT<Size>, SizeList*, SizesTCompare<Size>> SizesSizeMap;
typedef SizesSizeMap::iterator SizesSizeMapIter;

typedef map<SizesT<int>, SizeList*, SizesTCompare<int>> SizesIntMap;
typedef SizesIntMap::iterator SizesIntMapIter;

typedef vector<int> NumItersVec;

typedef map<SizesT<int>, const NumItersVec*, SizesTCompare<int>> SizesNumIterMap;
typedef SizesNumIterMap::iterator SizesNumIterMapIter;

struct RepSizesData {
  const SizeList *parent;
  const SizeList *controlParent;
  int size;
};

struct RepSizesCompare {
  bool operator() (const RepSizesData& lhs, const RepSizesData& rhs) const{
    if (lhs.parent < rhs.parent) {
      return true;
    }
    else if (lhs.parent == rhs.parent) {
      if (lhs.controlParent < rhs.controlParent) {
	return true;
      }
      else
	return lhs.size < rhs.size;
    }
    else
      return false;
  }
};

typedef map<RepSizesData, SizeList*, RepSizesCompare> RepSizesMap;
typedef RepSizesMap::iterator RepSizesMapIter;

typedef map<Size, SizeList*> SizeMap;
typedef SizeMap::iterator SizeMapIter;

typedef vector<const SizeList*> SizesVec;

struct RepeatedData {
  Size totSize;
  int reps;
};

struct RepeatedCompare {
  bool operator() (const RepeatedData& lhs, const RepeatedData& rhs) const{
    if (lhs.totSize < rhs.totSize) {
      return true;
    }
    else if (lhs.totSize == rhs.totSize) {
      return lhs.reps < rhs.reps;
    }
    else
      return false;
  }
};

typedef map<RepeatedData, SizeList*, RepeatedCompare> RepeatedMap;
typedef RepeatedMap::iterator RepeatedMapIter;


#if DOTENSORS
struct DistSizesData {
  const SizeList *parent;
  DistEntry entry;
};

struct DistSizesCompare {
  bool operator() (const DistSizesData& lhs, const DistSizesData& rhs) const{
    if (lhs.parent < rhs.parent) {
      return true;
    }
    else if (lhs.parent == rhs.parent) {
      return lhs.entry.m_val < rhs.entry.m_val;
    }
    else
      return false;
  }
};

typedef map<DistSizesData, SizeList*, DistSizesCompare> DistSizesMap;
typedef DistSizesMap::iterator DistSizesMapIter;
#endif

class SizesCache
{
 public:
  SizesSizeMap m_midSizesMap; //mid sizes
  SizesIntMap m_rangeMap; //first or last partition in range
  RepSizesMap m_repSizesMap; // repeated sizes based on other sizeList
  SizesVec m_otherSizes; // catchall sizes created elsewhere and given to cache
  DistSizesMap m_distSizesMap; // sizes with distribution coefficient
  SizesNumIterMap m_numItersMap; // num iterations map for range
  SizeMap m_constMap; // const sizes
  RepeatedMap m_constRepMap; // repeated sizes


#if DOTENSORS

#endif
  
  ~SizesCache();

  void TakeSize(SizeList *size);

  const SizeList* GetCachedMidSize(const SizeList *parent,
				   Size size);
  const SizeList* GetCachedRepeatedSize(const SizeList *parent,
					const SizeList *controlParent,
					const int size);
  const SizeList* GetCachedRepeatedSize(Size size,
					unsigned int numRepeats);					
  const SizeList* GetCachedRange(bool start,
			      const SizeList *parent,
				 int stride);
  const NumItersVec* GetNumItersVec(const SizeList *controlParent,
				    const int size);

  const SizeList* GetCachedDistSize(const SizeList *parent,
				    DistEntry entry);

  const SizeList* GetConstSize(Size size);
};
