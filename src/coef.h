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



#include <string>

using namespace std;

enum CoefVal { COEFVALALPHA,
	    COEFVALBETA,
	    COEFVALZERO,
	    COEFVALONE,
	    COEFVALONEHALF,
	       COEFVALNEGONE,
	    COEFVALNEGONEHALF };

class Coef
{
 public:
  CoefVal m_val;
 Coef(CoefVal val) : m_val(val) {}
  string BLISStr() const;
  string ElemStr() const;
  string TenStr() const;
  Coef operator*(const Coef &coef) const;
  Coef operator/(const Coef &coef) const;
  void operator=(const Coef &coef);
  bool operator==(const Coef &coef) const {return m_val == coef.m_val;}
  bool operator!=(const Coef &coef) const {return !(m_val == coef.m_val);}
};


extern Coef COEFALPHA;
extern Coef COEFBETA;
extern Coef COEFZERO;
extern Coef COEFONE;
extern Coef COEFONEHALF;
extern Coef COEFNEGONE;
extern Coef COEFNEGONEHALF;
