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

#include "avx.h"

#if DOLLDLA


string AVX::GlobalDeclarations() {
  
  string decls = "\n//--------------- AVX Declarations ----------------\n";
  decls += "typedef union { __m256i v; int i[8]; } intvec_reg;\n";
  decls += "typedef union { __m256i v; long int i[4]; } lintvec_reg;\n";
  decls += "\n";
  decls += "#define SLI 0xff00000000000000\n";
  decls += "#define SI 0xff00000\n";
  decls += "intvec_reg lsm00;\n";
  decls += "intvec_reg lsm01;\n";
  decls += "intvec_reg lsm02;\n";
  decls += "intvec_reg lsm03;\n";
  decls += "intvec_reg lsm04;\n";
  decls += "intvec_reg lsm05;\n";
  decls += "intvec_reg lsm06;\n";
  decls += "intvec_reg lsm07;\n";
  decls += "\n";
  decls += "lintvec_reg llsm00;\n";
  decls += "lintvec_reg llsm01;\n";
  decls += "lintvec_reg llsm02;\n";
  decls += "lintvec_reg llsm03;\n";

  decls += "//--------------- End AVX Declarations ----------------\n";
  return decls;
}

string AVX::SetupFunc() {
  string func = "void " + SetupFuncName() + "() {\n";
  for (int i = 0; i < 8; i++) {
    string ind = std::to_string((long long int) i);
    string rName = "lsm0" + ind;
    func += "\t" + rName + ".v = _mm256_setzero_si256();\n";
    for (int j = 0; j <= i; j++) {
      func += "\t" + rName + ".i[" + ind + "] = SI;\n";
    }
  }
  func += "\n";
  for (int i = 0; i < 4; i++) {
    string ind = std::to_string((long long int) i);
    string rName = "llsm0" + ind;
    func += "\t" + rName + ".v = _mm256_setzero_si256();\n";
    for (int j = 0; j <= i; j++) {
      func += "\t" + rName + ".i[" + ind + "] = SLI;\n";
    }
  }
  func += "}\n";
  return func;
}

#endif // DOLLDLA
