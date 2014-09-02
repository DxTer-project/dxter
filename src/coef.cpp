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



#include "coef.h"
#include <stdio.h>
#include <ostream>
#include <iostream>

Coef COEFALPHA(COEFVALALPHA);
Coef COEFBETA(COEFVALBETA);
Coef COEFZERO(COEFVALZERO);
Coef COEFONE(COEFVALONE);
Coef COEFTWO(COEFVALTWO);
Coef COEFONEHALF(COEFVALONEHALF);
Coef COEFNEGONE(COEFVALNEGONE);
Coef COEFNEGONEHALF(COEFVALNEGONEHALF);

string Coef::BLISStr() const
{
  switch(m_val)
    {
    case (COEFVALALPHA):
      return "alpha";
    case (COEFVALBETA):
      return "beta";
    case (COEFVALZERO):
      return "&BLIS_ZERO";
    case (COEFVALONE):
      return "&BLIS_ONE";
    case (COEFVALONEHALF):
      return "&BLIS_ONE_HALF";
    case (COEFVALNEGONE):
      return "&BLIS_MINUS_ONE";
    case (COEFVALNEGONEHALF):
      return "&BLIS_MINUS_ONE_HALF";
    default:
      throw;
    }
}

string Coef::ElemStr() const
{
  switch(m_val)
    {
    case (COEFVALALPHA):
      return "alpha";
    case (COEFVALBETA):
      return "beta";
    case (COEFVALZERO):
      return "0.0";
    case (COEFVALONE):
      return "1.0";
    case (COEFVALONEHALF):
      return "0.5";
    case (COEFVALNEGONE):
      return "-1.0";
    case (COEFVALNEGONEHALF):
      return "-0.5";
    default:
      throw;
    }
}

string Coef::TenStr() const
{
  switch(m_val)
    {
    case (COEFVALALPHA):
      return "alpha";
    case (COEFVALBETA):
      return "beta";
    case (COEFVALZERO):
      return "0.0";
    case (COEFVALONE):
      return "1.0";
    case (COEFVALTWO):
      return "2.0";
    case (COEFVALONEHALF):
      return "0.5";
    case (COEFVALNEGONE):
      return "-1.0";
    case (COEFVALNEGONEHALF):
      return "-0.5";
    default:
      throw;
    }
}

string Coef::LLDLAStr() const
{
  switch(m_val)
    {
    case (COEFVALALPHA):
      return "alpha";
    case(COEFVALBETA):
      return "beta";
    case(COEFVALZERO):
      return "&LLDLA_ZERO";
    case(COEFVALONE):
      return "&LLDLA_ONE";
    case(COEFVALONEHALF):
      return "&LLDLA_ONE_HALF";
    case(COEFVALNEGONE):
      return "&LLDLA_NEG_ONE";
    case(COEFVALNEGONEHALF):
      return "&LLDLA_NEG_ONE_HALF";
    default:
      throw;
    }
}

Coef Coef::operator*(const Coef &coef) const
{
  if (m_val == COEFVALONE)
    return coef;
  if (m_val == COEFVALZERO)
    return *this;
  if (coef.m_val == COEFVALONE)
    return *this;
  if (coef.m_val == COEFVALZERO)
    return coef;
  else {
    cout << "Coef*Coef = " << ElemStr() << " * " << coef.ElemStr() << endl;
    throw;
  }

}

Coef Coef::operator/(const Coef &coef) const
{
  if (m_val == COEFVALZERO)
    return *this;
  if (coef.m_val == COEFVALONE)
    return *this;
  if (coef.m_val == COEFVALZERO) {
    cout << "operator/ divide by zero\n";
    throw;
  }
  else {
    cout << "Coef/Coef = " << ElemStr() << " / " << coef.ElemStr() << endl;
    throw;
  }
}


void Coef::operator=(const Coef &coef)
{
  m_val = coef.m_val;
}
