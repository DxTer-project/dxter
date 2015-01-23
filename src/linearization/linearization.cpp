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
#include "clearLinElem.h"

Linearization::~Linearization()
{
  for(auto clear : m_clears) {
    delete clear;
  }
  m_clears.clear();
  m_order.clear();
  m_cost = 0;
}

void Linearization::InsertVecClearing()
{
  /*
    go through the list and find the diff between
    input and output variables.  add clears
    after for inputs that aren't output and aren't used later in the list

    ommit if a variable is still found in an outer linearization
  */
  throw;  
}

Cost Linearization::Cost() 
{
  if (m_cost <= 0) {
    throw;
  }
  else
    return m_cost;
}

void Linearization::operator=(const Linearization &rhs)
{
  if (m_clears.empty()) {
    m_order = rhs.m_order;
  }
  else {
    m_order.reserve(rhs.m_order.size());
    m_clears.reserve(rhs.m_clears.size());
    for(auto elem : rhs.m_order) {
      if (elem->IsClear()) {
	ClearLinElem *clear = new ClearLinElem(*((ClearLinElem*)elem));
	m_order.push_back(clear);
	m_clears.push_back(clear);
      }
      else {
	m_order.push_back(elem);
      }
    }
  }
}
