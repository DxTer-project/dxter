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



#include "layers.h"
#include "rqoTuple.h"
#include "rqoFieldValuePair.h"
#include "rqoBasis.h"
#include <unordered_map>
#include <iostream>

#if DORQO

bool satisfiesJoin(Tuple tuple1, Tuple tuple2, int key1, int key2)
{
    if(tuple1.getFields().at(key1).getValue() == 
        tuple2.getFields().at(key2).getValue())
    {
        return true;
    }
    return false;
}

Tuple joinTuples(Tuple one, Tuple two, int key)
{
    Tuple newTuple;
    for(FieldValuePair fvPair : one.getFields())
    {
        newTuple.addField(fvPair);
    }
    for(FieldValuePair fvPair : two.getFields())
    {
        if(fvPair.getValue() != two.getFields().at(key).getValue())
        {
            newTuple.addField(fvPair);
        }
    }
    return newTuple;
}


vector<Tuple> nestedJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2)
{
    vector<Tuple> output;


    for(auto tup1 : list1)
    {
        for(auto tup2 : list2)
        {
            if(satisfiesJoin(tup1, tup2, key1, key2))
            {
                Tuple toAdd = joinTuples(tup1, tup2, key2);
                output.push_back(toAdd);
            }
        }
    }

    return output;
}

vector<Tuple> hashJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2)
{
    unordered_map<string, Tuple> hash;
    vector<Tuple> output;
    unordered_map<string, Tuple>::const_iterator iter = hash.begin();

    for(auto tup1 : list1)
    {
        hash.insert(pair<string, Tuple>(
            tup1.getFields().at(key1).getValue(), tup1));
    }
    for(auto tup2 : list2)
    {
        iter = hash.find(tup2.getFields().at(key2).getValue());
        if(iter != hash.end())
        {
            Tuple toAdd = joinTuples(iter->second, tup2, key2);
            output.push_back(toAdd);
        }
    }
    return output;
}

vector<Tuple> mergeJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2)
{
    vector<Tuple> output;

    return output;
}


vector<Tuple> mergeFunc(vector<Tuple> left, vector<Tuple> right)
{
    vector<Tuple> output;

    while(left.size() != 0 && right.size() != 0)
    {
        if((*left.begin()).compareTo((*right.begin()), 0))
        {
            output.push_back((*left.begin()));
            left.erase(left.begin());
        }
        else
        {
            output.push_back((*right.begin()));
            right.erase(right.begin());
        }
    }
    while(left.size() != 0)
    {
        output.push_back((*left.begin()));
        left.erase(left.begin());
    }
    while(right.size() != 0)
    {
        output.push_back((*right.begin()));
        right.erase(right.begin());
    }

    return output;
}

vector<Tuple> sortFunc(vector<Tuple> list)
{
    if(list.size() <= 1)
    {
        return list;
    }

    int mid = list.size() / 2;
    vector<Tuple>::iterator iter = list.begin() + mid;
    vector<Tuple> left(list.begin(), iter);
    vector<Tuple> right(++iter, list.end());

    left = sortFunc(left);
    right = sortFunc(right);

    return mergeFunc(left, right);

}

vector<Tuple> unionFunc(vector<Tuple> list1, vector<Tuple> list2)
{
    list1 = sortFunc(list1);
    list2 = sortFunc(list2);

    vector<Tuple> output;

    while(list1.size() != 0 && list2.size() != 0)
    {
        if((*list1.begin()).equals((*list2.begin())))
        {
            output.push_back((*list1.begin()));
            list1.erase(list1.begin());
            list2.erase(list2.begin());
        }
        else if((*list1.begin()).compareTo((*list2.begin()), 0))
        {
            output.push_back((*list1.begin()));
            list1.erase(list1.begin());
        }
        else
        {
            output.push_back((*list2.begin()));
            list2.erase(list2.begin());
        }
    }
    while(list1.size() != 0)
    {
        output.push_back((*list1.begin()));
        list1.erase(list1.begin());
    }
    while(list2.size() != 0)
    {
        output.push_back((*list2.begin()));
        list2.erase(list2.begin());
    }

    return output;
}


#endif