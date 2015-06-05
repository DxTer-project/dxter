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

#include <ostream>
#include "coef.h"

using namespace std;

class Coef;

enum IndStreamType
  { BLISSTREAM,
    ELEMSTREAM,
    TENSORSTREAM,
    LLDLASTREAM,
    BOOLSTREAM,
    RQOSTREAM,
    OTHERSTREAM };

// Class to output code with indentations
class IndStream
{
 public:
  unsigned int m_level;
  ostream *o;
  IndStreamType m_type;
 IndStream(ostream *out, IndStreamType type) : m_level(0), o(out), m_type(type) {}
  void operator++();
  void operator--();
  void operator<<(const Coef &coef);
  void flush();
  ostream& operator*();
  void Indent(unsigned int offset=0);  
  string Tabs(unsigned int offset);
  string LoopLevel(unsigned int offset=0) const;
};
