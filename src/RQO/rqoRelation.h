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



#pragma once

#include "layers.h"
#include "rqoBasis.h"
#include "rqoAttribute.h"
#include "rqoTuple.h"
#include <unordered_map>

#if DORQO



class Relation
{
public:
    vector<Attribute> attributes;
    vector<Tuple> tuples;
    string m_name;

    Relation(string name) {m_name = name;}
    virtual void addAttribute(string name, string type, bool indexable);
    virtual void addTuple(Tuple tuple) {tuples.push_back(tuple);}
    virtual void printTable();
    virtual vector<Tuple> getTuples() {return tuples;}
    virtual string getName() {return m_name;}
    virtual int getSize() {return tuples.size();}
    virtual double getSelectivity(int key);
};


#endif
