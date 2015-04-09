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

#ifndef ISA_EXTENSION_H_
#define ISA_EXTENSION_H_

#include "transform.h"

#if DOLLDLA

class ISAExtension {
protected:
string m_name;
vector<pair<string, SingleTrans*>> m_archTrans;

public:
string SetupFuncName() { return m_name + "_setup"; }
vector<pair<string, SingleTrans*>> GetArchTrans() { return m_archTrans; }
  
virtual string SetupFunc() = 0;
virtual string GlobalDeclarations() = 0;
};

#endif // DOLLDLA

#endif
