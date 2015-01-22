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

#include "base.h"
#include "layers.h"

#if DOLLDLA

class Architecture
{
 public:
  // Single precision
  virtual int SVecRegWidth() = 0;
  virtual string SVecRegTypeDec() = 0;
  virtual string STypeName() = 0;
  virtual string SAddCode(string operand1, string operand2, string result) = 0;
  virtual string SMulCode(string operand1, string operand2, string result) = 0;
  virtual string SFMACode(string operand1, string operand2, string operand3, string result) = 0;
  virtual string SAccumCode(string memPtr, string startinLoc) = 0;
  virtual string SContiguousLoad(string memPtr, string receivingLoc) = 0;
  virtual string SStridedLoad(string memPtr, string receivingLoc, string stride) = 0;
  virtual string SDuplicateLoad(string memPtr, string receivingLoc) = 0;
  virtual string SContiguousStore(string memPtr, string startingLoc) = 0;
  virtual string SStridedStore(string memPtr, string startingLoc, string stride) = 0;
  virtual string SZeroVar(string varName) = 0;
  virtual double SFlopsPerCycle() = 0;

  // Double precision
  virtual int DVecRegWidth() = 0;
  virtual string DVecRegTypeDec() = 0;
  virtual string DTypeName() = 0;
  virtual string DAddCode(string operand1, string operand2, string result) = 0;
  virtual string DMulCode(string operand1, string operand2, string result) = 0;
  virtual string DFMACode(string operand1, string operand2, string operand3, string result) = 0;
  virtual string DAccumCode(string memPtr, string startinLoc) = 0;
  virtual string DContiguousLoad(string memPtr, string receivingLoc) = 0;
  virtual string DStridedLoad(string memPtr, string receivingLoc, string stride) = 0;
  virtual string DDuplicateLoad(string memPtr, string receivingLoc) = 0;
  virtual string DContiguousStore(string memPtr, string startingLoc) = 0;
  virtual string DStridedStore(string memPtr, string startingLoc, string stride) = 0;
  virtual string DZeroVar(string varName) = 0;
  virtual double DFlopsPerCycle() = 0;

  // Compilation
  virtual string CompileString(string executableName, string testFileName) = 0;

  // Performance
  virtual double CyclesPerSecond() = 0;

  // General
  int VecRegWidth(Type type);
  string TypeName(Type type);
  string VecRegTypeDec(Type type);
  string AddCode(Type type, string operand1, string operand2, string result);
  string MulCode(Type type, string operand1, string operand2, string result);
  string FMACode(Type type, string operand1, string operand2, string operand3, string result);
  string AccumCode(Type type, string memPtr, string startingLoc);
  string ContiguousLoad(Type type, string memPtr, string receivingLoc);
  string StridedLoad(Type type, string memPtr, string receivingLoc, string stride);
  string DuplicateLoad(Type type, string memPtr, string receivingLoc);
  string ContiguousStore(Type type, string memPtr, string startingLoc);
  string StridedStore(Type type, string memPtr, string startingLoc, string stride);
  string ZeroVar(Type type, string varName);
  double FlopsPerCycle(Type type);

  // Cost model
  virtual int ContigVecLoadCost() = 0;
  virtual int ContigVecStoreCost() = 0;
};

class AMDEngSample : public Architecture
{
 public:
  // Single precision
  virtual int SVecRegWidth();
  virtual string SVecRegTypeDec();
  virtual string STypeName();
  virtual string SAddCode(string operand1, string operand2, string result);
  virtual string SMulCode(string operand1, string operand2, string result);
  virtual string SFMACode(string operand1, string operand2, string operand3, string result);
  virtual string SAccumCode(string memPtr, string startingLoc);
  virtual string SContiguousLoad(string memPtr, string receivingLoc);
  virtual string SStridedLoad(string memPtr, string receivingLoc, string stride);
  virtual string SDuplicateLoad(string memPtr, string receivingLoc);
  virtual string SContiguousStore(string memPtr, string startingLoc);
  virtual string SStridedStore(string memPtr, string startingLoc, string stride);
  virtual string SZeroVar(string varName);
  virtual double SFlopsPerCycle();

  // Double precision
  virtual int DVecRegWidth();
  virtual string DVecRegTypeDec();
  virtual string DTypeName();
  virtual string DAddCode(string operand1, string operand2, string result);
  virtual string DMulCode(string operand1, string operand2, string result);
  virtual string DFMACode(string operand1, string operand2, string operand3, string result);
  virtual string DAccumCode(string memPtr, string startinLoc);
  virtual string DContiguousLoad(string memPtr, string receivingLoc);
  virtual string DStridedLoad(string memPtr, string receivingLoc, string stride);
  virtual string DDuplicateLoad(string memPtr, string receivingLoc);
  virtual string DContiguousStore(string memPtr, string startingLoc);
  virtual string DStridedStore(string memPtr, string startingLoc, string stride);
  virtual string DZeroVar(string varName);
  virtual double DFlopsPerCycle();

  virtual string CompileString(string executableName, string testFileName);
  virtual double CyclesPerSecond();

  // Cost model
  virtual int ContigVecLoadCost();
  virtual int ContigVecStoreCost();
};

class Stampede : public Architecture
{
 public:
  // Single precision
  virtual int SVecRegWidth();
  virtual string SVecRegTypeDec();
  virtual string STypeName();
  virtual string SAddCode(string operand1, string operand2, string result);
  virtual string SMulCode(string operand1, string operand2, string result);
  virtual string SFMACode(string operand1, string operand2, string operand3, string result);
  virtual string SAccumCode(string memPtr, string startingLoc);
  virtual string SContiguousLoad(string memPtr, string receivingLoc);
  virtual string SStridedLoad(string memPtr, string receivingLoc, string stride);
  virtual string SDuplicateLoad(string memPtr, string receivingLoc);
  virtual string SContiguousStore(string memPtr, string startingLoc);
  virtual string SStridedStore(string memPtr, string startingLoc, string stride);
  virtual string SZeroVar(string varName);
  virtual double SFlopsPerCycle();

  // Double precision
  virtual int DVecRegWidth();
  virtual string DVecRegTypeDec();
  virtual string DTypeName();
  virtual string DAddCode(string operand1, string operand2, string result);
  virtual string DMulCode(string operand1, string operand2, string result);
  virtual string DFMACode(string operand1, string operand2, string operand3, string result);
  virtual string DAccumCode(string memPtr, string startinLoc);
  virtual string DContiguousLoad(string memPtr, string receivingLoc);
  virtual string DStridedLoad(string memPtr, string receivingLoc, string stride);
  virtual string DDuplicateLoad(string memPtr, string receivingLoc);
  virtual string DContiguousStore(string memPtr, string startingLoc);
  virtual string DStridedStore(string memPtr, string startingLoc, string stride);
  virtual string DZeroVar(string varName);
  virtual double DFlopsPerCycle();

  // Compilation
  virtual string CompileString(string executableName, string testFileName);
  virtual double CyclesPerSecond();

  // Cost model
  virtual int ContigVecLoadCost();
  virtual int ContigVecStoreCost();
};

class HaswellMacbook : public Stampede
{
 public:
  // Single precision
  virtual string SFMACode(string operand1, string operand2, string operand3, string result);
  virtual double SFlopsPerCycle();

  // Double precision
  virtual string DFMACode(string operand1, string operand2, string operand3, string result);
  virtual double DFlopsPerCycle();

  // Compilation
  virtual string CompileString(string executableName, string testFileName);
  virtual double CyclesPerSecond();

  // Cost model
  virtual int ContigVecLoadCost();
  virtual int ContigVecStoreCost();
};

extern Architecture* arch;

#endif // DOLLDLA
