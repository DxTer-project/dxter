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



#include "IndStream.h"
#include <sstream>

void IndStream::operator++()
{
  ++m_level;
}

void IndStream::operator--()
{
  --m_level;
}

ostream& IndStream::operator*()
{
  return *o;
}

void IndStream::flush()
{
  o->flush();
}

void IndStream::Indent(unsigned int offset)
{
  for(unsigned int i = 0; i < (offset+m_level); ++i) {
    *o << "\t";
  }
}

string IndStream::Tabs(unsigned int offset)
{
  string ret;
  for(unsigned int i = 0; i < (offset+m_level); ++i) {
    ret += "\t";
  }
  return ret;
}

string IndStream::LoopLevel(unsigned int offset) const
{
  std::stringstream str;
  str << m_level + offset;
  return str.str();
}

void IndStream::operator<<(const Coef &coef)
{
  if (m_type == BLISSTREAM) {
    *(*this) << coef.BLISStr();
  }
  else if (m_type == ELEMSTREAM) {
    *(*this) << coef.ElemStr();
  }
  else {
    *(*this) << coef.TenStr();
  }
}
