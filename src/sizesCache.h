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

using namespace std;

template <typename T>
struct SizesT {
  const Sizes *parent;
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


typedef map<SizesT<Size>, const Sizes*, SizesTCompare<Size>> SizesSizeMap;
typedef SizesSizeMap::iterator SizesSizeMapIter;

typedef map<SizesT<int>, const Sizes*, SizesTCompare<int>> SizesIntMap;
typedef SizesIntMap::iterator SizesIntMapIter;

typedef vector<int> NumItersVec;

typedef map<SizesT<int>, const NumItersVec*, SizesTCompare<int>> SizesNumIterMap;
typedef SizesNumIterMap::iterator SizesNumIterMapIter;

struct RepSizesData {
  const Sizes *parent;
  const Sizes *controlParent;
  Size size;
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

typedef map<RepSizesData, const Sizes*, RepSizesCompare> RepSizesMap;
typedef RepSizesMap::iterator RepSizesMapIter;




class SizesCache
{
 public:
  SizesSizeMap m_midSizesMap;
  SizesIntMap m_rangeStartMap;
  SizesIntMap m_rangeEndMap;
  RepSizesMap m_repSizesMap;
  SizesNumIterMap m_numItersMap;

  ~SizesCache();

  const Sizes* GetCachedMidSize(const Sizes *parent,
				Size size,
				const NumItersVec *numIters);
  const Sizes* GetCachedRepeatedSize(const Sizes *parent,
				     const Sizes *controlParent,
				     const Size size,
				     const NumItersVec *numIters);
  const Sizes* GetCachedRange(bool start,
			      const Sizes *parent,
			      int stride,
			      const NumItersVec *numIters);
  const NumItersVec* GetNumItersVec(const Sizes *controlParent,
				    const Size size);
};
