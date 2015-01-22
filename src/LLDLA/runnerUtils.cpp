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

#include "runnerUtils.h"

#if DOLLDLA

unique_ptr<ImplementationMap> ImpStrMap(Universe* uni) {
  std::unique_ptr<ImplementationMap> impMap(new ImplementationMap());
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

string NoWhitespace(const string str) {
  for (auto ch : str) {
    if (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r') {
      cout << "ERROR: Unallowed whitespace in string: " << str;
      throw;
    }
  }
  return str;
}

string UtilsIncludesAndDefines() {
  return "#include <math.h>\n#include <pmmintrin.h>\n#include <stdio.h>\n#include <stdlib.h>\n#include <time.h>\n#include <string.h>\n#define VEC_PD_LOAD(ptr) _mm_load_pd((ptr))\n#define VEC_DUP_LOAD(ptr) _mm_loaddup_pd((ptr))\n#define VEC_D_LOAD(ptr1, ptr2) _mm_loadh_pd(_mm_load_sd((ptr1)), (ptr2));\n#define VEC_PD_STORE(ptr, vec_reg) _mm_store_pd((ptr), (vec_reg))\ntypedef union	{\n__m128d v;\ndouble d[2];\n} v2df_t;\n";
}

string Rdtsc() {
  return "unsigned long long rdtsc() {unsigned long long int x;unsigned a, d;__asm__ volatile(\"rdtsc\" : \"=a\" (a), \"=d\" (d));return ((unsigned long long)a) | (((unsigned long long)d) << 32);}";
}

string UtilFuncs() {
  cout << "ERROR: UtilFuncs() is not implemented yet." << endl;
  throw;
  return "";
}

string Utils() {
  string utils = UtilsIncludesAndDefines();
  utils += Rdtsc();
  utils += UtilFuncs();
  return utils;
}

#endif // DOLLDLA
