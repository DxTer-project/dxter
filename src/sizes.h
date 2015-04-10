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



#include <vector>
#include <math.h>
#include "IndStream.h"

using namespace std;

typedef double Size;
typedef double Cost;


enum SizesType { MIDSIZES,
		 REPEATEDSIZES,
		 RANGESIZES,
		 BADSIZE };
	       

/*
  for each type:
  MIDSIZES -> m_valA is val
           -> m_valB is totalSize
  REPEATEDSIZES -> m_valA is val
                -> m_valC is num repeats
  RANGESIZES -> m_valA is start
             -> m_valB is end
             -> m_valC is stride  
*/

//Iterator for Sizes
class SizesIter
{
 private:
  Size Update(Size size) const;
 public:
  int m_currPos;
  Size m_valA, m_valB;
  int m_valC;
  SizesType m_type;
  double m_coeff;
  int m_repeatNum;
  int m_repeats;
 SizesIter() : m_currPos(-1) {}
  SizesIter(Size valA, Size valB, int valC, SizesType type, double coeff, int repeats);
  void operator++();
  Size operator*() const;
  void operator=(const SizesIter &rhs);
  bool AtEnd() const;
  bool AtCurrRepeatEnd() const;
};

class SizeEntry
{
 public:
  Size m_valA, m_valB;
  int m_valC;
  SizesType m_type;
  int m_repeats;
  SizeEntry();

  SizesType GetType() const;
  bool ConstantMidSizeEntry() const;



  void SetRepeatedSizes(Size size, int repeats);
  void SetSizeRange(Size start, int stride, Size end);
  void SetMidSizes(Size size, Size totalSize);
  Size operator[] (unsigned int n) const;
  void operator= (const SizeEntry &rhs);
  bool operator==(const SizeEntry &rhs) const;
  bool operator!=(const SizeEntry &rhs) const;
  bool operator==(const Size &rhs) const;
  bool operator!=(const Size &rhs) const;
  bool operator<= (const Size &rhs) const;
  bool EvenlyDivisibleBy(const Size &size) const;
  void Print() const;
  void Print(IndStream &out) const;
  unsigned int NumSizesPerRepeat() const;
  unsigned int NumSizes() const;
  bool IsZero() const;
  SizesIter GetIter(double coeff) const;
  Cost Sum() const;
  Cost SumSquares() const;
  SizeEntry SumWith(const SizeEntry &rhs) const;
};

typedef vector<SizeEntry*> EntryVec;
typedef EntryVec::iterator EntryVecIter;
typedef EntryVec::const_iterator EntryVecConstIter;

//This class keeps track of the sizes of all
// matrices that flow on a wire
class Sizes
{
 private:
  Size Update(Size size) const;
 public:

  double m_constVal;
  double m_coeff;
  EntryVec m_entries;
  Sizes();
  Sizes(double constVal);
  Sizes(const Sizes &rhs);
  virtual ~Sizes();

  Sizes* PartitionedSize(int splitPoint);

  void Print() const;
  void Print(IndStream &out) const;
  void AddRepeatedSizes(Size size, int repeats);
  void AddSizesWithLimit(Size start, int stride, Size end);
  void AddMidSizes(Size size, Size totalSize);
  void ClearSizes();
  void SetCoeff(double coeff);
  unsigned int NumSizes() const;
  Size operator[] (unsigned int n) const;
  void operator= (const Sizes &rhs);
  bool operator== (const Sizes &rhs) const;
  bool operator!= (const Sizes &rhs) const;
  bool operator== (const Size &rhs) const;
  bool operator!= (const Size &rhs) const;
  bool operator<= (const Size &rhs) const;
  bool operator>(const Size &rhs) const;
  bool EvenlyDivisibleBy(const Size &size) const;
  Cost Sum() const;
  Cost SumSquares() const;
  Cost SumCubes() const;
  Cost SumProds11(const Sizes &sizes) const;
  Cost SumProds21(const Sizes &sizes) const;
  Cost SumProds111(const Sizes &sizes1, const Sizes &sizes2) const;
  void PairwiseSum(const Sizes &sizes1, const Sizes &sizes2);
  bool AllOnes() const;
  SizesIter GetIter(unsigned int sizeNum) const;
  bool IsZero(unsigned int n) const;

  bool IsPartitionable(const Size partitionPoint) const;
  bool IsConstant() const;
  Size OnlyEntry() const;
};

extern Sizes *ONES;

typedef Sizes *SizesArray;
