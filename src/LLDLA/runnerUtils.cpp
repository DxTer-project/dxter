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


string TypeToStr(Type type) {
  switch(type) {
  case(REAL_SINGLE):
    return "real_single_precision";
  case(REAL_DOUBLE):
    return "real_double_precision";
  default:
    cout << "ERROR: Bad type in TypeToStr" << endl;
    LOG_FAIL("replacement for throw call");
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
    LOG_FAIL("replacement for throw call");
  }
}

string NoWhitespace(const string str) {
  for (auto ch : str) {
    if (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r') {
      cout << "ERROR: Unallowed whitespace in string: " << str;
      LOG_FAIL("replacement for throw call");
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
  LOG_FAIL("replacement for throw call");
  return "";
}

string Utils() {
  string utils = UtilsIncludesAndDefines();
  utils += Rdtsc();
  utils += UtilFuncs();
  return utils;
}

#endif // DOLLDLA
