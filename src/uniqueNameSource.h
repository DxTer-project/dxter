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

#ifndef UNIQUE_NAME_SOURCE_H_
#define UNIQUE_NAME_SOURCE_H_

#include <string>

#include "uniqueNumSource.h"

using namespace std;

#if DOLLDLA

class UniqueNameSource {
 private:
  UniqueNumSource m_numSource;
  string m_prefix;

 public:
  UniqueNameSource(string prefix);
  string Next(string name);
};

extern UniqueNameSource* localInputNames;

#endif // DOLLDLA

#endif
