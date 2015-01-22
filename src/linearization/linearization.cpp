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

#include "linearization.h"

Linearization::~Linearization()
{
  for(auto clear : m_clears) {
    delete clear;
  }
  m_clears.clear();
  m_order.clear();
}

void InsertVecClearing()
{
  /*
    go through the list and find the diff between
    input and output variables.  add clears
    after for inputs that aren't output and aren't used later in the list

    ommit if a variable is still found in an outer linearization
  */
  throw;  
}
