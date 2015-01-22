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



#include "TrProps.h"
#include "string.h"
#include "base.h"

TrProps::TrProps(bool invert, Side side, Tri tri, Diag diag, Trans trans, Coef coeff, Type type)
  :m_invert(invert), m_side(side), m_tri(tri), m_diag(diag), m_trans(trans), m_coeff(coeff), m_type(type)
{

}

void TrProps::Duplicate(const TrProps *trmm)
{
  m_invert = trmm->m_invert;
  m_side = trmm->m_side;
  m_tri = trmm->m_tri;
  m_diag = trmm->m_diag;
  m_trans = trmm->m_trans;
  m_coeff = trmm->m_coeff;
  m_type = trmm->m_type;
}

void TrProps::FlattenCore(ofstream &out) const
{
  WRITE(m_invert);
  WRITE(m_side);
  WRITE(m_tri);
  WRITE(m_diag);
  WRITE(m_trans);
  WRITE(m_coeff);
  WRITE(m_type);
}


void TrProps::UnflattenCore(ifstream &in, SaveInfo &info) 
{
  READ(m_invert);
  READ(m_side);
  READ(m_tri);
  READ(m_diag);
  READ(m_trans);
  READ(m_coeff);
  READ(m_type);
}


