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

#include "logging.h"
#include "sizes.h"
#include "sizesCache.h"

using namespace std;

#define COMBINESIZEENTRIES 1

SizeList ONEVAL(1);
SizeList *ONES = &ONEVAL;

SizesCache SizeList::M_cache;

SizesIter::SizesIter(Size valA, Size valB, int valC, SizesType type, double coeff, int repeats)
:m_currPos(0),
m_valA(valA),
m_valB(valB),
m_valC(valC),
m_type(type),
 m_coeff(coeff),
 m_repeatNum(0),
 m_repeats(repeats)
{  
}

void SizesIter::operator++()
{
  if (m_currPos < 0) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  ++m_currPos;
  if (AtCurrRepeatEnd()) {
    ++m_repeatNum;
    m_currPos = 0;
  }
}

Size SizesIter::operator*() const
{
  if (m_currPos < 0) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (AtEnd()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (m_type == REPEATEDSIZES) {
    return Update(m_valA);
  }
  else if (m_type == MIDSIZES) {
    int numFullIters = (int)floor(m_valB / m_valA);
    if (numFullIters > m_currPos)
      return Update(m_valA);
    else if (numFullIters+1 > m_currPos)
      return Update(m_valB - numFullIters*m_valA);
    else {
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
  else if (m_type == RANGESIZES) {
    return Update(m_valA + m_currPos * m_valC);
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
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
  m_repeatNum = rhs.m_repeatNum;
  m_repeats = rhs.m_repeats;
}

bool SizesIter::AtEnd() const
{
  return ((m_repeatNum == m_repeats) ||
      (m_repeatNum == m_repeats-1 &&
       AtCurrRepeatEnd()));
}

bool SizesIter::AtCurrRepeatEnd() const
{
  if (m_type == MIDSIZES) {
    if (m_valB == 0) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    if (m_currPos == 0) {
      return false;
    }
    double size = (double)m_valB;
    double iters = ceil(size / (double)m_valA);
    return m_currPos >= iters;
  }
  else if (m_type == REPEATEDSIZES) {
    return m_currPos >= ((double)m_valC);
  }
  else if (m_type == RANGESIZES) {
    //Add one so if it's something like 0:B:n*B we calc n+1 iterations not n
    //the corner case of it being 0:B:n*B-1 is unlikely, so this is faster
    double sizeDiff = abs((double)m_valB - m_valA);
    double eachPortion = sizeDiff +1;
    double iters = ceil(eachPortion / abs((double)m_valC));
    return m_currPos >= iters;
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

SizeEntry::SizeEntry()
  : m_type(BADSIZE), m_repeats(-1)
{
}

SizesType SizeEntry::GetType() const {
  return m_type;
}

bool SizeEntry::ConstantMidSizeEntry() const {
  return (((int) m_valB) % ((int) m_valA)) == 0;
}

void SizeEntry::SetRepeatedSizes(Size size, int repeats)
{
  m_valA = size;
  m_valC = repeats;
  m_type = REPEATEDSIZES;
  m_repeats = 1;
}

void SizeEntry::SetSizeRange(Size start, int stride, Size end)
{
  if (start == 0 && end == 0) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  m_valA = start;
  m_valB = end;
  m_valC = stride;
  m_type = RANGESIZES;
  m_repeats = 1;
}

void SizeEntry::SetMidSizes(Size size, Size totalSize)
{
  //  cout << "ERROR: calling SetMidSizes" << endl;
  //  LOG_FAIL("replacement for throw call");
  m_valA = size;
  m_valB = totalSize;
  m_type = MIDSIZES;
  m_repeats = 1;
}

void SizeEntry::Print() const
{
  cout << m_repeats << " x ";
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
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

void SizeEntry::Print(IndStream &out) const
{
  *out << m_repeats << " x ";
  if (m_type == REPEATEDSIZES)
    *out << "RepeatedSize " << m_valA
	 << ", " << m_valC << " times\n";
  else if (m_type == MIDSIZES)
    *out << "MidSize " << m_valA
    << ", total of " << m_valB
    << endl;
  else if (m_type == RANGESIZES)
    *out << "RangeSize " << m_valA << ":"
    << m_valC << ":"
    << m_valB << endl;
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}
  

bool SizeEntry::operator==(const SizeEntry &rhs) const
{
  if (m_repeats != rhs.m_repeats)
    return false;
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
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

Size SizeEntry::operator[] (unsigned int n) const
{
  if ((int)(n / NumSizesPerRepeat()) > m_repeats) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  n = n % NumSizesPerRepeat();
  unsigned int iter = n;
  if (m_type == REPEATEDSIZES) {
    if ((int)iter >= (double)m_valC) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    return m_valA;
  }
  else if (m_type == MIDSIZES) {
    double adjustedValB;
    
    adjustedValB = (double)m_valB;
    
    unsigned int numFullIters = (int)floor(adjustedValB / m_valA);
    if (numFullIters > iter)
      return m_valA;
    else if (numFullIters+1 > iter)
      return adjustedValB - numFullIters*m_valA;
    else {
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
  else if (m_type == RANGESIZES) {
    if (iter >= NumSizes()) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    return m_valA + ((Size)iter) * m_valC;
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

void SizeEntry::operator= (const SizeEntry &rhs)
{
  m_valA = rhs.m_valA;
  m_valB = rhs.m_valB;
  m_valC = rhs.m_valC;
  m_type = rhs.m_type;
  m_repeats = rhs.m_repeats;
}

bool SizeEntry::operator!=(const SizeEntry &rhs) const
{
  return !(*this == rhs);
}

bool SizeEntry::operator==(const Size &rhs) const
{
  switch(m_type)
    {
    case (MIDSIZES):
      {
	if (m_valA != rhs)
	  return false;
	return (fmod(m_valB, m_valA) == 0);
      }
    case (REPEATEDSIZES):
      {
	return m_valA == rhs;
      }
    case (RANGESIZES):
      {
	return false;
	break;
      }
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

bool SizeEntry::EvenlyDivisibleBy(const Size &size) const
{
  switch(m_type)
    {
    case (MIDSIZES):
      {
	if (!fmod(m_valA,size))
	  return false;
	return (fmod(m_valB, m_valA) == 0);
      }
    case (REPEATEDSIZES):
      {
	return !fmod(m_valA,size);
      }
    case (RANGESIZES):
      {
	return !fmod(m_valA, size) 
	  && !fmod(m_valC, size);
	break;
      }
    default:
      LOG_FAIL("replacement for throw call");
      throw;
    }
}

bool SizeEntry::operator!=(const Size &rhs) const
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
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

unsigned int SizeEntry::NumSizesPerRepeat() const
{
  if (m_type == MIDSIZES) {
    return (unsigned int)ceil(m_valB / (double)m_valA);
  }
  else if (m_type == REPEATEDSIZES) {
    return (unsigned int)m_valC;
  }
  else if (m_type == RANGESIZES) {
    //Add one so if it's something like 0:B:n*B we calc n+1 iterations not n
    double sizeDiff = abs((double)m_valB - m_valA);
    double eachPortion = sizeDiff+1;
    double num = ceil(eachPortion / abs((double)m_valC));
    return (unsigned int)num;
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
}

unsigned int SizeEntry::NumSizes() const
{
  return NumSizesPerRepeat() * m_repeats;
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
  SizesIter iter(m_valA,m_valB,m_valC,m_type,coeff,m_repeats);
  return iter;
}

SizeEntry SizeEntry::SumWith(const SizeEntry &rhs) const
{
  if (NumSizes() != rhs.NumSizes()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (m_repeats != rhs.m_repeats) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
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
    else {
      LOG_FAIL("replacement for throw call");
      throw;
    }
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
    else {
      return rhs.SumWith(*this);

    }
  }
  else if (m_type == RANGESIZES) {
    if (rhs.m_type == RANGESIZES) {
      entry.SetSizeRange(rhs.m_valA+m_valA,
                         rhs.m_valC+m_valC,
                         rhs.m_valB+m_valB);
      return entry;
    }
    else {
      return rhs.SumWith(*this);
    }
  }
  else {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  LOG_FAIL("replacement for throw call");
  throw;
}

SizeList::SizeList()
{
  m_constVal = NAN;
  m_coeff = 0;
  m_cached = false;
}

SizeList::SizeList(double constVal)
{
  m_constVal = constVal;
  m_coeff = 0;
  m_cached = false;
}

SizeList::SizeList(const SizeList &rhs)
{
  m_constVal = NAN;
  m_coeff = 0;
  *this = rhs;
  m_cached = false;
}


SizeList::~SizeList()
{
  if (IsCached())
    throw;
  EntryVecIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter)
    delete *iter;
  m_entries.clear();
}


void SizeList::Print() const
{
  cout << m_entries.size() << " size entries\n";
  cout << NumSizes() << " sizes\n";
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    cout << "With coef " << m_coeff << ", ";
    (*iter)->Print();
  }
}

void SizeList::Print(IndStream &out) const
{
  *out << m_entries.size() << " size entries\n";
  *out << NumSizes() << " sizes\n";
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    (*iter)->Print(out);
  }
}

void SizeList::AddRepeatedSizes(Size size, int repeats)
{
  if (IsCached())
    throw;
#if COMBINESIZEENTRIES
  if (!m_entries.empty()) {
    SizeEntry *oldEntry = m_entries.back();
    if (oldEntry->m_type == REPEATEDSIZES) {
      if (oldEntry->m_valA == size) {
	if (oldEntry->m_valC == repeats) {
	  oldEntry->m_repeats++;
	  return;
	}
      }
    }
  }
#endif
  SizeEntry *entry = new SizeEntry;
  entry->SetRepeatedSizes(size, repeats);
  m_entries.push_back(entry);
}

unsigned int SizeList::AddMidSizes(Size size, Size totalSize)
{
  if (IsCached())
    throw;
#if COMBINESIZEENTRIES
  if (!fmod(totalSize, size)){
    AddRepeatedSizes(size, totalSize / size);
    return totalSize / size;
  }
  if (!m_entries.empty()) {
    SizeEntry *oldEntry = m_entries.back();
    if (oldEntry->m_type == MIDSIZES) {
      if (oldEntry->m_valA == size) {
	if (oldEntry->m_valB == totalSize) {
	  oldEntry->m_repeats++;
	  return oldEntry->NumSizesPerRepeat();
	}
      }
    }
  }
#endif
  SizeEntry *entry = new SizeEntry;
  entry->SetMidSizes(size, totalSize);
  m_entries.push_back(entry);
  return entry->NumSizesPerRepeat();
}

unsigned int SizeList::AddSizesWithLimit(Size start, int stride, Size end)
{
  if (IsCached())
    throw;
  if (stride == 0) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
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
  if (start == 0 && end == 0) {
    AddRepeatedSizes(0, 1);
    return 1;
  }
  else {
#if COMBINESIZEENTRIES
    if (!m_entries.empty()) {
      SizeEntry *oldEntry = m_entries.back();
      if (oldEntry->m_type == RANGESIZES) {
	if (oldEntry->m_valA == start) {
	  if (oldEntry->m_valB == end) {
	    if (oldEntry->m_valC == stride) {
	      oldEntry->m_repeats++;
	      return oldEntry->NumSizesPerRepeat();
	    }
	  }
	}
      }
    }
#endif
    SizeEntry *entry = new SizeEntry;
    entry->SetSizeRange(start, stride, end);
    m_entries.push_back(entry);
    return entry->NumSizesPerRepeat();
  }
}

void SizeList::SetCoeff(double coeff)
{
  if (IsCached())
    throw;
  if (coeff == 0) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (coeff == 1)
    return;
  m_coeff = coeff;
}

void SizeList::ClearSizes()
{
  if (IsCached())
    throw;
  EntryVecIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter)
    delete *iter;
  m_entries.clear();
}

Size SizeList::operator[] (unsigned int n) const
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
  cout << n << " of " << NumSizes() << endl;
  LOG_FAIL("replacement for throw call");
  throw;
}


void SizeList::operator=(const SizeList &rhs)
{
  if (IsCached())
    throw;

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

bool SizeList::operator==(const SizeList &rhs) const
{
  if ((!std::isnan((double)m_constVal) && std::isnan((double)(rhs.m_constVal)))
      || (std::isnan((double)m_constVal) && !std::isnan((double)(rhs.m_constVal)))) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (!std::isnan((double)m_constVal) && (m_constVal != rhs.m_constVal))
    return false;
  if (m_entries.size() != rhs.m_entries.size()) {
    cout << m_entries.size() << endl;
    cout << rhs.m_entries.size() << endl;
    Print();
    cout << "****\n";
    rhs.Print();
    LOG_FAIL("replacement for throw call");
    throw;
  }
  EntryVecConstIter iter1 = m_entries.begin();
  EntryVecConstIter iter2 = rhs.m_entries.begin();
  for(; iter1 != m_entries.end(); ++iter1,++iter2) {
    if (**iter1 != **iter2)
      return false;
  }
  return true;
}


bool SizeList::EffectivelyEqual(const SizeList &rhs) const
{
  if ((!std::isnan((double)m_constVal) && std::isnan((double)(rhs.m_constVal)))
      || (std::isnan((double)m_constVal) && !std::isnan((double)(rhs.m_constVal)))) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (!std::isnan((double)m_constVal) && (m_constVal != rhs.m_constVal))
    return false;
  if (m_entries.size() == rhs.m_entries.size()) {
    if (*this == rhs)
      return true;
  }
  if (NumSizes() != rhs.NumSizes()) {
    cout << m_entries.size() << endl;
    cout << rhs.m_entries.size() << endl;
    Print();
    cout << "****\n";
    rhs.Print();
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (!NumSizes())
    return true;

  EntryVecConstIter iter1 = m_entries.begin();
  EntryVecConstIter iter2 = rhs.m_entries.begin();
  SizesIter iterLeft = (*iter1)->GetIter(m_coeff);
  SizesIter iterRight = (*iter2)->GetIter(rhs.m_coeff);
  bool keepGoing = true;
  while (keepGoing) {
    if (iterLeft.AtEnd()) {
      ++iter1;
      if (iter1 == m_entries.end()) {
	break;
      }
      iterLeft = (*iter1)->GetIter(m_coeff);
    }
    if (iterRight.AtEnd()) {
      ++iter2;
      if (iter2 == m_entries.end()) {
	break;
      }
      iterRight = (*iter2)->GetIter(rhs.m_coeff);
    }
    if (*iterLeft != *iterRight)
      return false;
    ++iterLeft;
    ++iterRight;
  }
  return true;
}

bool SizeList::operator!=(const SizeList &rhs) const
{
  return !(*this == rhs);
}

bool SizeList::operator==(const Size &rhs) const
{
  if (!std::isnan((double)m_constVal) && (m_constVal != rhs))
    return false;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    if (**iter != rhs)
      return false;
  }
  return true;
}

bool SizeList::EvenlyDivisibleBy(const Size &size) const
{
  if (!std::isnan((double)m_constVal))
    return !(fmod(m_constVal,size));
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    if (!(*iter)->EvenlyDivisibleBy(size))
      return false;
  }
  return true;

}


bool SizeList::operator!=(const Size &rhs) const
{
  return !(*this == rhs);
}

unsigned int SizeList::NumSizes() const
{
  unsigned int num = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    num += (*iter)->NumSizes();
  }
  return num;
}

Cost SizeList::Sum() const
{
  if (!std::isnan((double)m_constVal)) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  Cost cost = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    const SizeEntry *entry = *iter;
    if (entry->m_type == REPEATEDSIZES) {
      cost += entry->m_repeats * entry->m_valC
	* Update(entry->m_valA);
    }
    else if (entry->m_type == MIDSIZES) {
      double numIters = entry->m_valB/((double)(entry->m_valA));
      int numFullIters = (int)floor(numIters);
      cost += entry->m_repeats * numFullIters * Update(entry->m_valA);
    }
    else if (entry->m_type == RANGESIZES) {
      SizesIter iter2 = entry->GetIter(m_coeff);
      while(!iter2.AtEnd()) {
        cost += *iter2;
        ++iter2;
      }
    }
    else {
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
  return cost;
}

bool SizeList::IsZero(unsigned int n) const
{
  return m_entries[n]->IsZero();
}

Cost SizeList::SumSquares() const
{
  if (!std::isnan((double)m_constVal)) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  Cost cost = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    const SizeEntry *entry = *iter;
    if (entry->m_type == REPEATEDSIZES) {
      cost += entry->m_repeats * entry->m_valC * pow(Update(entry->m_valA),2);
    }
    else if (entry->m_type == MIDSIZES) {
      double numIters = entry->m_valB/((double)(entry->m_valA));
      int numFullIters = (int)floor(numIters);
      cost += entry->m_repeats * numFullIters * pow(Update(entry->m_valA),2);
    }
    else if (entry->m_type == RANGESIZES) {
      SizesIter iter2 = entry->GetIter(m_coeff);
      while(!iter2.AtEnd()) {
        cost += pow(*iter2,2);
        ++iter;
      }
    }
    else {
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
  return cost;
}

Size SizeList::Update(Size size) const
{
  if (m_coeff)
    return ceil(size * m_coeff);
  else
    return size;
}

Cost SizeList::SumCubes() const
{
  if (!std::isnan((double)m_constVal)) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  Cost cost = 0;
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    const SizeEntry *entry = *iter;
    if (entry->m_type == REPEATEDSIZES) {
      cost += entry->m_repeats * (double)entry->m_valC
	* pow(Update(entry->m_valA),3);
    }
    else if (entry->m_type == MIDSIZES) {
      double numIters = entry->m_valB/((double)(entry->m_valA));
      int numFullIters = (int)floor(numIters);
      cost += entry->m_repeats * numFullIters
	* pow(Update(entry->m_valA),3);
    }
    else if (entry->m_type == RANGESIZES) {
      SizesIter iter2 = entry->GetIter(m_coeff);
      while(!iter2.AtEnd()) {
        cost += pow(*iter2,3);
        ++iter;
      }
    }
    else {
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
  return cost;
}

Cost SizeList::SumProds11(const SizeList &sizes) const
{
  if (!std::isnan((double)m_constVal)) {
    return m_constVal * sizes.Sum();
  }
  if (!std::isnan((double)(sizes.m_constVal))) {
    return sizes.m_constVal * Sum();
  }
  if (NumSizes() != sizes.NumSizes()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  Cost cost = 0;
  unsigned int j = 0;
  SizesIter iter2 = sizes.GetIter(j);
  for(unsigned int i = 0; i < m_entries.size(); ++i) {
    SizesIter iter1 = GetIter(i);
    while (!iter1.AtEnd()) {
      if (iter2.AtEnd()) {
	++j;
	iter2 = sizes.GetIter(j);
      }
      cost += *iter1 * *iter2;
      ++iter1;
      ++iter2;
    }
  }
  return cost;
}

Cost SizeList::SumProds21(const SizeList &sizes) const
{
  if (!std::isnan((double)m_constVal)) {
    return m_constVal * m_constVal * sizes.Sum();
  }
  if (!std::isnan((double)(sizes.m_constVal))) {
    return sizes.m_constVal * SumSquares();
  }
  if (m_entries.size() != sizes.m_entries.size()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (NumSizes() != sizes.NumSizes()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  Cost cost = 0;
  for(unsigned int i = 0; i < m_entries.size(); ++i) {
    SizesIter iter1 = GetIter(i);
    SizesIter iter2 = sizes.GetIter(i);
    while (!iter1.AtEnd()) {
      cost += *iter1 * *iter1 * *iter2;
      ++iter1;
      ++iter2;
    }
    if (!iter2.AtEnd()) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
  return cost;
}

Cost SizeList::SumProds111(const SizeList &sizes1, const SizeList &sizes2) const
{
  if (!std::isnan((double)m_constVal)) {
    return m_constVal * sizes1.SumProds11(sizes2);
  }
  if (!std::isnan((double)(sizes1.m_constVal))) {
    return sizes1.m_constVal * SumProds11(sizes2);
  }
  if (!std::isnan((double)(sizes2.m_constVal))) {
    return sizes2.m_constVal * SumProds11(sizes1);
  }
  if (m_entries.size() != sizes1.m_entries.size() || m_entries.size() != sizes2.m_entries.size()) {
    if (NumSizes() != sizes1.NumSizes() || NumSizes() != sizes2.NumSizes()) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
    else {
      Cost cost = 0;
      for (unsigned int i = 0; i < NumSizes(); ++i) {
	cost += (*this)[i] * sizes1[i] * sizes2[i];
      }
      return cost;
    }
  }
  Cost cost = 0;
  for(unsigned int i = 0; i < m_entries.size(); ++i) {
    SizesIter iter1 = GetIter(i);
    SizesIter iter2 = sizes1.GetIter(i);
    SizesIter iter3 = sizes2.GetIter(i);
    while (!iter1.AtEnd()) {
      Size size1 = *iter1;
      Size size2 = *iter2;
      Size size3 = *iter3;
      cost += size1 * size2 * size3;
      ++iter1;
      ++iter2;
      ++iter3;
    }
    if (!iter2.AtEnd()) {
      cout << iter1.m_currPos << endl;
      cout << iter2.m_currPos << endl;
      cout << iter3.m_currPos << endl;
      LOG_FAIL("replacement for throw call");
      throw;
    }
    if (!iter3.AtEnd()) {
      LOG_FAIL("replacement for throw call");
      throw;
    }
  }
  return cost;
}

bool SizeList::AllOnes() const
{
  if (m_constVal == 1)
    return true;
  if (!m_entries.size()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
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

SizesIter SizeList::GetIter(unsigned int sizeNum) const
{
  if (sizeNum >= m_entries.size()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  
  return m_entries[sizeNum]->GetIter(m_coeff);
}

void SizeList::PairwiseSum(const SizeList &sizes1, const SizeList &sizes2)
{
  if (sizes1.NumSizes() != sizes2.NumSizes()) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (sizes1.m_coeff || sizes2.m_coeff) {
    LOG_FAIL("replacement for throw call");
    throw;
  }
  if (!std::isnan((double)(sizes1.m_constVal))) {
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
          int numIters = (int)ceil(entry->m_valB/((double)(entry->m_valA)));
          newEntry->SetMidSizes(entry->m_valA+constVal,entry->m_valB+constVal*numIters);
	  newEntry->m_repeats = entry->m_repeats;
        }
        else if (entry->m_type == REPEATEDSIZES) {
          newEntry->SetRepeatedSizes(entry->m_valA+constVal,entry->m_valC);
	  newEntry->m_repeats = entry->m_repeats;
        }
        else if (entry->m_type == RANGESIZES) {
          newEntry->SetSizeRange(entry->m_valA+constVal,
                                 entry->m_valC+constVal,
                                 entry->m_valB+constVal);
	  newEntry->m_repeats = entry->m_repeats;
        }
        else {
          LOG_FAIL("replacement for throw call");
	  throw;
	}
        m_entries.push_back(newEntry);
      }
    }
  }
  if (!std::isnan((double)(sizes2.m_constVal))) {
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

bool SizeList::operator<= (const Size &rhs) const
{
  EntryVecConstIter iter = m_entries.begin();
  for(; iter != m_entries.end(); ++iter) {
    if (! (*(*iter) <= rhs))
      return false;
  }
  return true;
}

bool SizeList::operator>(const Size &rhs) const
{
  return !(*this <= rhs);
}

bool SizeList::IsConstant() const {

  for (auto entry : this->m_entries) {
    SizesType entryType = entry->GetType();
    if (entryType == RANGESIZES) {
      cout << "FOUND RANGESIZES" << endl;
      return false;
    } else if (entryType == MIDSIZES) {
      if (!(entry->ConstantMidSizeEntry())) {
	return false;
      }
    }
  }
  return true;
}

Size SizeList::OnlyEntry() const {
  if (!IsConstant()) {
    cout << "ERROR: sizes not a constant" << endl;
    cout << "# entries = " << m_entries.size()  << endl;
    LOG_FAIL("replacement for throw call");
    throw;
  }
  return (*this)[0];
}

const SizeList* GetConst(Size val)
{
  return SizeList::M_cache.GetConstSize(val);
}
