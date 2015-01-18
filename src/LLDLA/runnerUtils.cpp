/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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

#include "runnerUtils.h"

#if DOLLDLA

unique_ptr<ImplementationMap> ImpStrMap(Universe* uni) {
  std::unique_ptr<ImplementationMap> impMap = std::make_unique<ImplementationMap>();
  GraphNum i;
  for (i = 1; i <= uni->TotalCount(); i++) {
    std::stringbuf sbuf;
    std::ostream out(&sbuf);
    IndStream istream = IndStream(&out, LLDLASTREAM);
    uni->Print(istream, i);
    impMap->insert(NumImplementationPair(i, sbuf.str()));
  }
  return impMap;
}

string TypeToStr(Type type) {
  switch(type) {
  case(REAL_SINGLE):
    return "real_single_precision";
  case(REAL_DOUBLE):
    return "real_double_precision";
  default:
    cout << "ERROR: Bad type in TypeToStr" << endl;
    throw;
  }
}

string VecTypeToStr(VecType vecType) {
  switch(vecType) {
  case(ROWVECTOR):
    return "row_vector";
  case(COLVECTOR):
    return "column_vector";
  default:
    cout << "ERROR: Bad type in TypeToStr" << endl;
    throw;
  }
}

string DateAndTimeString() {
  time_t t = time(0);
  struct tm* now = localtime(&t);
  return std::to_string((long long int) (now->tm_year + 1900)) + "-"
    + std::to_string(now->tm_mon + 1) + "-"
    + std::to_string(now->tm_mday) + "-"
    + std::to_string(now->tm_hour) + "-"
    + std::to_string(now->tm_min) + "-"
    + std::to_string(now->tm_sec);
}

string NoWhitespace(string str) {
  for (auto ch : str) {
    if (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r') {
      cout << "ERROR: Unallowed whitespace in string: " << str;
      throw;
    }
  }
  return str;
}

#endif // DOLLDLA
