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

#pragma once

#include "DLANode.h"
#include "DLAOp.h"
#include "LLDLA.h"

#if DOLLDLA

class FMAdd : public DLAOp<3, 1>
{
public:
 virtual NodeType GetType() const { return "FMAdd"; }
 static Node* BlankInst() { return new FMAdd; }
 virtual Node* GetNewInst() { return BlankInst(); }
 
 virtual void Prop();
 virtual void PrintCode(IndStream &out);
 virtual ClassType GetNodeClass() const { return GetClass(); }
 static ClassType GetClass() { return "FMAdd"; }

 virtual bool IsReadOnly() const { return false; }
 virtual bool IsDataDependencyOfInput() const { return true; }
};

class Add : public DLAOp<2, 1>
{
public:
 virtual NodeType GetType() const { return "Add"; }
 static Node* BlankInst() { return new Add; }
 virtual Node* GetNewInst() { return BlankInst(); }
 
 virtual void Prop();
 virtual void PrintCode(IndStream &out);
 virtual ClassType GetNodeClass() const { return GetClass(); }
 static ClassType GetClass() { return "Add"; }

 virtual bool IsReadOnly() const { return false; }
 virtual bool IsDataDependencyOfInput() const { return true; }
};

class Mul : public DLAOp<2, 1>
{
public:
 virtual NodeType GetType() const { return "Mul"; }
 static Node* BlankInst() { return new Mul; }
 virtual Node* GetNewInst() { return BlankInst(); }
 
 virtual void Prop();
 virtual void PrintCode(IndStream &out);
 virtual ClassType GetNodeClass() const { return GetClass(); }
 static ClassType GetClass() { return "Mul"; }

 virtual bool IsReadOnly() const { return false; }
 virtual bool IsDataDependencyOfInput() const { return true; }
};

#endif // DOLLDLA
