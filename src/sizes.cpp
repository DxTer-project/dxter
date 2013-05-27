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



#include "sizes.h"
#include <cmath>
#include <cstring>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <ostream>
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <map>
#include <list>
#include <fstream>
#include <cstdlib>

using namespace std;

Sizes ONEVAL(1);
Sizes *ONES = &ONEVAL;

SizesIter::SizesIter(Size valA, Size valB, int valC, SizesType type, double coeff)
:m_currPos(0),
m_valA(valA),
m_valB(valB),
m_valC(valC),
m_type(type),
m_coeff(coeff)
{
}

void SizesIter::operator++()
{
  if (m_currPos < 0)
    throw;
  ++m_currPos;
}

Size SizesIter::operator*() const
{
  if (m_currPos < 0)
    throw;
  if (m_type == REPEATEDSIZES) {
    if (m_currPos >= m_valC)
      throw;
    return Update(m_valA);
  }
  else if (m_type == MIDSIZES) {
    int numFullIters = floor(m_valB / m_valA);
    if (numFullIters > m_currPos)
      return Update(m_valA);
    else if (numFullIters+1 > m_currPos)
      return Update(m_valB - numFullIters*m_valA);
    else
      throw;
  }
  else if (m_type == RANGESIZES) {
    if (m_currPos >= ceil((double)abs(m_valA-m_valB)))
      throw;
    return Update(m_valA + m_currPos * m_valC);
  }
  else
    throw;
}

Size SizesIter::Update(Size size) const
{
  if (m_coeff)
    return ceil(m_coeff * size);
  else
    return size;
}
void SizesIter::operator=(const SizesIter &rhs)
{
  m_valA = rhs.m_valA;
  m_valB = rhs.m_valB;
  m_valC = rhs.m_valC;
  m_currPos = rhs.m_currPos;
  m_coeff = rhs.m_coeff;
  m_type = rhs.m_type;
}

bool SizesIter::AtEnd() const
{
  if (m_type == MIDSIZES) {
    if (m_valB == 0)
      throw;
    if (m_currPos == 0)
      return false;
    double done = m_currPos * m_valA;
    return done >= m_valB;
  }
  else if (m_type == REPEATEDSIZES) {
    return (m_currPos >= m_valC);
  }
  else if (m_type == RANGESIZES) {
    Size val = m_valA + m_valC * m_currPos;
    if (m_valC < 0) {
      return val < m_valB;
    }
    else
      return val > m_valB;
  }
  else
    throw;
}

SizeEntry::SizeEntry()
: m_type(BADSIZE)
{
}

void SizeEntry::SetRepeatedSizes(Size size, int repeats)
{
  m_valA = size;
  m_valC = repeats;
  m_type = REPEATEDSIZES;
}

void SizeEntry::SetSizeRange(Size start, int stride, Size end)
{
  if (start == 0 && end == 0)
    throw;
  m_valA = start;
  m_valB = end;
  m_valC = stride;
  m_type = RANGESIZES;
}

void SizeEntry::SetMidSizes(Size size, Size totalSize)
{
  m_valA = size;
  m_valB = totalSize;
  m_type = MIDSIZES;
}

void SizeEntry::Print() const
{
  if (m_type == REPEATEDSIZES)
    cout << "RepeatedSize " << m_valA
    << ", " << m_valC << " times\n";
  else if (m_type == MIDSIZES)
    cout << "MidSize " << m_valA
    << ", total of " << m_valB
    << endl;
  else if (m_type == RANGESIZES)
    cout << "RangeSize " << m_valA << ":"
    << m_valC << ":"
    << m_valB << endl;
  else
    throw;
}

bool SizeEntry::operator==(const SizeEntry &rhs) const
{
  if (m_type != rhs.m_type)
    return false;
  if (m_type == MIDSIZES) {
    return m_valA == rhs.m_valA &&
    m_valB == rhs.m_valB;
  }
  else if (m_type == REPEATEDSIZES) {
    return m_valA == rhs.m_valA &&
    m_valC == rhs.m_valC;
    
  }
  else if (m_type == RANGESIZES) {
    return m_valA == rhs.m_valA &&
    m_valB == rhs.m_valB &&
    m_valC == rhs.m_valC;
  }
  else
    throw;
}

Size SizeEntry::operator[] (unsigned int n) const
{
  if (m_type == REPEATEDSIZES) {
    if ((int)n >= m_valC)
      throw;
    return m_valA;
  }
  else if (m_type == MIDSIZES) {
    unsigned int numFullIters = floor(m_valB / m_valA);
    if (numFullIters > n)
      return m_valA;
    else if (numFullIters+1 > n)
      return m_valB - numFullIters*m_valA;
    else
      throw;
  }
  else if (m_type == RANGESIZES) {
    if (n >= NumSizes())
      throw;
    return m_valA + ((Size)n) * m_valC;
  }
  else
    throw;
}

void SizeEntry::operator= (const SizeEntry &rhs)
{
  m_valA = rhs.m_valA;
  m_valB = rhs.m_valB;
  m_valC = rhs.m_valC;
  m_type = rhs.m_type;
}

bool SizeEntry::operator!=(const SizeEntry &rhs) const
{
  return !(*this == rhs);
}

bool SizeEntry::operator<= (const Size &rhs) const
{
  switch(m_type) {
  case (MIDSIZES):
    return m_valA <= rhs;
  case (REPEATEDSIZES):
    return m_valA <= rhs;
  case (RANGESIZES):
    {
      if (m_valC > 0)
	return m_valB <= rhs;
      else
	return m_valA <= rhs;
    }
  default:
    throw;
  }
}

unsigned int SizeEntry::NumSizes() const
{
  if (m_type == MIDSIZES) {
    return ceil(((double)m_valB) / m_valA);
  }
  else if (m_type == REPEATEDSIZES) {
    return m_valC;
  }
  else if (m_type == RANGESIZES) {
    double diff = abs(m_valA - m_valB);
    return ceil(diff / abs(m_valC))+1;
  }
  else
    throw;
}

bool SizeEntry::IsZero() const
{
  if (m_type == REPEATEDSIZES)
    return m_valA == 0;
  else
    return false;
}

SizesIter SizeEntry::GetIter(double coeff) const
{
  SizesIter iter(m_valA,m_valB,m_valC,m_type,coeff);
  return iter;
}

SizeEntry SizeEntry::SumWith(const SizeEntry &rhs) const
{
  if (NumSizes() != rhs.NumSizes())
    throw;
  SizeEntry entry;
  if (m_type == MIDSIZES) {
    if (rhs.m_type == MIDSIZES) {
      entry.SetMidSizes(m_valA+rhs.m_valA,
                        m_valB+rhs.m_valB);
      return entry;
    }
    else if (rhs.m_type == REPEATEDSIZES) {
      entry.SetMidSizes(m_valA+rhs.m_valA,
                        m_valB+rhs.m_valA*rhs.m_valC);
      return entry;
    }
    else if (rhs.m_type == RANGESIZES) {
      entry.SetSizeRange(rhs.m_valA+m_valA,
                         rhs.m_valC,
                         rhs.m_valB+m_valA);
      return entry;
    }
    else
      throw;
  }
  else if (m_type == REPEATEDSIZES) {
    if (rhs.m_type == REPEATEDSIZES) {
      entry.SetRepeatedSizes(m_valA+rhs.m_valA,
                             m_valC);
      return entry;
    }
    else if (rhs.m_type == RANGESIZES) {
      entry.SetSizeRange(rhs.m_valA+m_valA,
                         rhs.m_valC,
                         rhs.m_valB+m_valA);
      return entry;
    }
    else
      return rhs.SumWith(*this);
  }
  else if (m_type == RANGESIZES) {
    if (rhs.m_type == RANGESIZES) {
      entry.SetSizeRange(rhs.m_valA+m_valA,
                         rhs.m_valC+m_valC,
                         rhs.m_valB+m_valB);
      return entry;
    }
    else
      return rhs.SumWith(*this);
  }
  else
    throw;
  throw;
}

Sizes::Sizes()
{
  m_constVal = NAN;
  m_coeff = 0;
}

Sizes::Sizes(double constVal)
{
  m_constVal = constVal;
  m_coeff = 0;
}

Sizes::~Sizes()
{
  EntryVecIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter)
    delete *iter;
  m_entries.clear();
}


void Sizes::Print() const
{
  cout << m_entries.size() << " size entries\n";
  cout << NumSizes() << " sizes\n";
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    (*iter)->Print();
  }
}

void Sizes::AddRepeatedSizes(Size size, int repeats)
{
  SizeEntry *entry = new SizeEntry;
  entry->SetRepeatedSizes(size, repeats);
  m_entries.push_back(entry);
}

void Sizes::AddMidSizes(Size size, Size totalSize)
{
  SizeEntry *entry = new SizeEntry;
  entry->SetMidSizes(size, totalSize);
  m_entries.push_back(entry);
}

void Sizes::AddSizesWithLimit(Size start, int stride, Size end)
{
  if (stride == 0)
    throw;
  if (stride < 0) {
    int mod = (int)start % (-1*stride);
    if (mod)
      start -= mod;
    else
      start += stride;
  }
  else {
    int mod = (int)end % stride;
    if (mod)
      end -= mod;
    else
      end -= stride;
  }
  SizeEntry *entry = new SizeEntry;
  if (start == 0 && end == 0)
    entry->SetRepeatedSizes(0, 1);
  else
    entry->SetSizeRange(start, stride, end);
  m_entries.push_back(entry);
}

void Sizes::SetCoeff(double coeff)
{
  if (coeff == 0)
    throw;
  if (coeff == 1)
    return;
  m_coeff = coeff;
}

void Sizes::ClearSizes()
{
  EntryVecIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter)
    delete *iter;
  m_entries.clear();
}

Size Sizes::operator[] (unsigned int n) const
{
  unsigned int currLoc = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    int numSizes = (*iter)->NumSizes();
    if (numSizes + currLoc > n) {
      return Update((*(*iter))[n-currLoc]);
    }
    else
      currLoc += numSizes;
  }
  throw;
}


void Sizes::operator=(const Sizes &rhs)
{
  m_coeff = rhs.m_coeff;
  m_constVal = rhs.m_constVal;
  EntryVecIter iter2 = m_entries.begin();
  for(; iter2 != m_entries.end(); ++iter2)
    delete *iter2;
  m_entries.clear();
  m_entries.reserve(rhs.m_entries.size());
  EntryVecConstIter iter = rhs.m_entries.begin();
  for(; iter != rhs.m_entries.end(); ++iter) {
    SizeEntry *entry = new SizeEntry;
    *entry = **iter;
    m_entries.push_back(entry);
  }
}

bool Sizes::operator==(const Sizes &rhs) const
{
  if ((!isnan(m_constVal) && isnan(rhs.m_constVal))
      || (isnan(m_constVal) && !isnan(rhs.m_constVal)))
    throw;
  if (!isnan(m_constVal) && (m_constVal != rhs.m_constVal))
    return false;
  if (m_entries.size() != rhs.m_entries.size())
    throw;
  EntryVecConstIter iter1 = m_entries.begin();
  EntryVecConstIter iter2 = rhs.m_entries.begin();
  for(; iter1 != m_entries.end(); ++iter1,++iter2) {
    if (**iter1 != **iter2)
      return false;
  }
  return true;
}

bool Sizes::operator!=(const Sizes &rhs) const
{
  return !(*this == rhs);
}

unsigned int Sizes::NumSizes() const
{
  unsigned int num = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    num += (*iter)->NumSizes();
  }
  return num;
}


Cost Sizes::Sum() const
{
  if (!isnan(m_constVal))
    throw;
  Cost cost = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    const SizeEntry *entry = *iter;
    if (entry->m_type == REPEATEDSIZES) {
      cost += entry->m_valC * Update(entry->m_valA);
    }
    else if (entry->m_type == MIDSIZES) {
      double numIters = entry->m_valB/((double)(entry->m_valA));
      int numFullIters = floor(numIters);
      int numPartIters = ceil(numIters)-numFullIters;
      cost += numFullIters * Update(entry->m_valA);
      if (numPartIters)
        cost += Update(entry->m_valB - numFullIters * entry->m_valA);
    }
    else if (entry->m_type == RANGESIZES) {
      SizesIter iter2 = entry->GetIter(m_coeff);
      while(!iter2.AtEnd()) {
        cost += *iter2;
        ++iter2;
      }
    }
    else
      throw;
  }
  return cost;
}

bool Sizes::IsZero(unsigned int n) const
{
  return m_entries[n]->IsZero();
}

Cost Sizes::SumSquares() const
{
  if (!isnan(m_constVal))
    throw;
  Cost cost = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    const SizeEntry *entry = *iter;
    if (entry->m_type == REPEATEDSIZES) {
      cost += entry->m_valC * pow(Update(entry->m_valA),2);
    }
    else if (entry->m_type == MIDSIZES) {
      double numIters = entry->m_valB/((double)(entry->m_valA));
      int numFullIters = floor(numIters);
      int numPartIters = ceil(numIters)-numFullIters;
      cost += numFullIters * Update(entry->m_valA);
      if (numPartIters)
        cost += pow(Update(entry->m_valB - numFullIters * entry->m_valA),2);
    }
    else if (entry->m_type == RANGESIZES) {
      SizesIter iter2 = entry->GetIter(m_coeff);
      while(!iter2.AtEnd()) {
        cost += pow(*iter2,2);
        ++iter;
      }
    }
    else
      throw;
  }
  return cost;
}

Size Sizes::Update(Size size) const
{
  if (m_coeff)
    return ceil(size * m_coeff);
  else
    return size;
}

Cost Sizes::SumCubes() const
{
  if (!isnan(m_constVal))
    throw;
  Cost cost = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    const SizeEntry *entry = *iter;
    if (entry->m_type == REPEATEDSIZES) {
      cost += entry->m_valC * pow(Update(entry->m_valA),3);
    }
    else if (entry->m_type == MIDSIZES) {
      double numIters = entry->m_valB/((double)(entry->m_valA));
      int numFullIters = floor(numIters);
      int numPartIters = ceil(numIters)-numFullIters;
      cost += numFullIters * Update(entry->m_valA);
      if (numPartIters)
        cost += pow(Update(entry->m_valB - numFullIters * entry->m_valA),3);
    }
    else if (entry->m_type == RANGESIZES) {
      SizesIter iter2 = entry->GetIter(m_coeff);
      while(!iter2.AtEnd()) {
        cost += pow(*iter2,3);
        ++iter;
      }
    }
    else
      throw;
  }
  return cost;
}

Cost Sizes::SumProds11(const Sizes &sizes) const
{
  if (!isnan(m_constVal)) {
    return m_constVal * sizes.Sum();
  }
  if (!isnan(sizes.m_constVal)) {
    return sizes.m_constVal * Sum();
  }
  if (m_entries.size() != sizes.m_entries.size())
    throw;
  Cost cost = 0;
  for(unsigned int i = 0; i < m_entries.size(); ++i) {
    if (IsZero(i) || sizes.IsZero(i))
      continue;
    SizesIter iter1 = GetIter(i);
    SizesIter iter2 = sizes.GetIter(i);
    while (!iter1.AtEnd()) {
      if (iter2.AtEnd()) {
        cout << "bad\n";
        throw;
      }
      cost += *iter1 + *iter2;
      ++iter1;
      ++iter2;
    }
    if (!iter2.AtEnd())
      throw;
  }
  return cost;
}

Cost Sizes::SumProds21(const Sizes &sizes) const
{
  if (!isnan(m_constVal)) {
    return m_constVal * m_constVal * sizes.Sum();
  }
  if (!isnan(sizes.m_constVal)) {
    return sizes.m_constVal * SumSquares();
  }
  if (m_entries.size() != sizes.m_entries.size())
    throw;
  if (NumSizes() != sizes.NumSizes())
    throw;
  Cost cost = 0;
  for(unsigned int i = 0; i < m_entries.size(); ++i) {
    SizesIter iter1 = GetIter(i);
    SizesIter iter2 = sizes.GetIter(i);
    while (!iter1.AtEnd()) {
      cost += *iter1 * *iter1 * *iter2;
      ++iter1;
      ++iter2;
    }
    if (!iter2.AtEnd())
      throw;
  }
  return cost;
}

Cost Sizes::SumProds111(const Sizes &sizes1, const Sizes &sizes2) const
{
  if (!isnan(m_constVal)) {
    return m_constVal * sizes1.SumProds11(sizes2);
  }
  if (!isnan(sizes1.m_constVal)) {
    return sizes1.m_constVal * SumProds11(sizes2);
  }
  if (!isnan(sizes2.m_constVal)) {
    return sizes2.m_constVal * SumProds11(sizes1);
  }
  if (m_entries.size() != sizes1.m_entries.size() || m_entries.size() != sizes2.m_entries.size())
    throw;
  Cost cost = 0;
  for(unsigned int i = 0; i < m_entries.size(); ++i) {
    SizesIter iter1 = GetIter(i);
    SizesIter iter2 = sizes1.GetIter(i);
    SizesIter iter3 = sizes2.GetIter(i);
    while (!iter1.AtEnd()) {
      cost += *iter1 * *iter2 * *iter3;
      ++iter1;
      ++iter2;
      ++iter3;
    }
    if (!iter2.AtEnd()) {
      cout << iter1.m_currPos << endl;
      cout << iter2.m_currPos << endl;
      cout << iter3.m_currPos << endl;
      throw;
    }
    if (!iter3.AtEnd())
      throw;
  }
  return cost;
}

bool Sizes::AllOnes() const
{
  if (m_constVal == 1)
    return true;
  if (!m_entries.size())
    throw;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    if ((*iter)->m_type == REPEATEDSIZES) {
      if ((*iter)->m_valA != 1)
        return false;
    }
    else
      return false;
  }
  return true;
}

SizesIter Sizes::GetIter(unsigned int sizeNum) const
{
  if (sizeNum >= m_entries.size())
    throw;
  
  return m_entries[sizeNum]->GetIter(m_coeff);
}

void Sizes::PairwiseSum(const Sizes &sizes1, const Sizes &sizes2)
{
  if (sizes1.NumSizes() != sizes2.NumSizes())
    throw;
  if (sizes1.m_coeff || sizes2.m_coeff)
    throw;
  if (!isnan(sizes1.m_constVal)) {
    if (sizes1.m_constVal == 1) {
      *this = sizes2;
      return;
    }
    else {
      EntryVecConstIter iter = sizes2.m_entries.begin();
      for(; iter != sizes2.m_entries.end(); ++iter) {
        SizeEntry *entry = *iter;
        SizeEntry *newEntry = new SizeEntry;
        double constVal = sizes1.m_constVal;
        if (entry->m_type == MIDSIZES) {
          int numIters = ceil(entry->m_valB/((double)(entry->m_valA)));
          newEntry->SetMidSizes(entry->m_valA+constVal,entry->m_valB+constVal*numIters);
        }
        else if (entry->m_type == REPEATEDSIZES) {
          newEntry->SetRepeatedSizes(entry->m_valA+constVal,entry->m_valC);
        }
        else if (entry->m_type == RANGESIZES) {
          newEntry->SetSizeRange(entry->m_valA+constVal,
                                 entry->m_valC+constVal,
                                 entry->m_valB+constVal);
        }
        else
          throw;
        m_entries.push_back(newEntry);
      }
    }
  }
  if (!isnan(sizes2.m_constVal)) {
    PairwiseSum(sizes2, sizes1);
    return;
  }
  EntryVecConstIter iter1 = sizes1.m_entries.begin();
  EntryVecConstIter iter2 = sizes2.m_entries.begin();
  for(; iter1 != sizes1.m_entries.end(); ++iter1, ++iter2) {
    SizeEntry *entry = new SizeEntry;
    *entry = (*iter1)->SumWith(**iter2);
    m_entries.push_back(entry);
  }
}

bool Sizes::operator<= (const Size &rhs) const
{
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    if (! (*(*iter) <= rhs))
      return false;
  }
  return true;
}
