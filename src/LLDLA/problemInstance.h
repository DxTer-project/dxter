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

#ifndef PROBLEM_INSTANCE_H_
#define PROBLEM_INSTANCE_H_

#include "LLDLA.h"

#if DOLLDLA

class ProblemInstance {
 private:
  string* m_name;
  Type m_type;
  Cost m_cost;

 public:
  ProblemInstance();
  ~ProblemInstance();

  void AddDimension(int val, string dimName);

  Cost GetCost();
  string GetName();
  Type GetType();

  void SetCost(Cost cost);
  void SetName(string name);
  void SetType(Type type);
};

#endif // DOLLDLA

#endif // PROBLEM_INSTANCE_H_
