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

#include "problemInstance.h"

#if DOLLDLA

ProblemInstance::ProblemInstance() {
  m_name = new string("");
}

ProblemInstance::~ProblemInstance() {
  delete m_name;
}

void ProblemInstance::AddDimension(int val, string dimName) {
  string* oldName = m_name;
  m_name = new string(*oldName + "_" + dimName + std::to_string((long long int) val));
  delete oldName;
}

Cost ProblemInstance::GetCost() {
  return m_cost;
}

string ProblemInstance::GetName() {
  return *m_name;
}

Type ProblemInstance::GetType() {
  return m_type;
}

void ProblemInstance::SetCost(Cost cost) {
  m_cost = cost;
}

void ProblemInstance::SetName(string name) {
  delete m_name;
  m_name = new string(name);
}

void ProblemInstance::SetType(Type type) {
  m_type = type;
}

#endif // DOLLDLA
