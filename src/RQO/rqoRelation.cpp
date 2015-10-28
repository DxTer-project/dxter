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



#include "rqoRelation.h"

#if DORQO

void Relation::printTable()
{
    cout << m_name << " contains the following attributes:" << endl;
    for(auto attri : attributes)
    {
        cout << attri.m_name << " of type " << attri.m_type << endl;
    }
    for(auto tuple : tuples)
    {
        tuple.printTuple();
    }
}

double Relation::getSelectivity(int key)
{
    unordered_map<string, Tuple> hash;
    unordered_map<string, Tuple>::const_iterator iter = hash.begin();

    for(auto tup : tuples)
    {
        iter = hash.find(tup.getValueAt(key));
        if(iter == hash.end())
        {
            hash.insert(pair<string, Tuple>(
            tup.getValueAt(key), tup));
        }
    }

    return 1/hash.size();
}

#endif
