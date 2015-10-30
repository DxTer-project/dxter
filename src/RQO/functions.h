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
#include "rqoTuple.h"
#include "rqoFieldValuePair.h"
#include "rqoBasis.h"
#include "queryNodes.h"
#include "rqoRelation.h"
#include <unordered_map>
#include <iostream>
#include <algorithm>

using namespace queryNodes;

#if DORQO

vector<Tuple> scanFunc(Relation table, FieldValue query);
vector<Tuple> indexFunc(Relation table, OrNode query, int index);
vector<Tuple> nindexFunc(Relation table, OrNode query, set<int> indeces);
vector<Tuple> orderedindexFunc(Relation table, OrNode query, int index);
bool satisfiesJoin(Tuple tuple1, Tuple tuple2, int key1, int key2);
Tuple joinTuples(Tuple one, Tuple two, int key);
vector<Tuple> nestedJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2);
vector<Tuple> hashJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2);
vector<Tuple> mergeJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2);
vector<Tuple> leftOuterJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2);
vector<Tuple> rightOuterJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2);
vector<Tuple> mergeFunc(vector<Tuple> left, vector<Tuple> right, int key1, int key2);
vector<Tuple> sortFunc(vector<Tuple> list, int key);
vector<Tuple> unionFunc(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2);
vector<Tuple> crossProduct(vector<Tuple> list1, vector<Tuple> list2);
vector<Tuple> projection(vector<Tuple> list, vector<string> values);
vector<Tuple> fullOuterJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2);
void printTuples(vector<Tuple> list);

#endif