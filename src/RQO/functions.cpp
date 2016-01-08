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



#include "functions.h"

#if DORQO

vector<Tuple> scanFunc(Relation table, OrNode query, vector<string> values)
{
	vector<Tuple> output;
    cout << "in scan" << endl;

	for(auto tuple : table.getTuples())
	{
        cout << "gonna evaluate" << endl;
		if(query.evaluate(tuple, -1))
		{
            cout << "evaluate success?" << endl;
			output.push_back(tuple);
		}
	}

    //set up newvalues
    vector<string> newValues;
    for(auto fvPair : output.at(0).getFields())
    {
        newValues.push_back(fvPair.getField());
    }
    vector<string>::iterator it = newValues.begin();
	for(auto string : values)
    {
        it = find(newValues.begin(), newValues.end(), string);
        if(it != newValues.end())
        {
            newValues.erase(it);
        }
        
    }

    output = trim(output, newValues);

	return output;
}

vector<Tuple> indexFunc(Relation table, OrNode query, int index, vector<string> values)
{
    vector<Tuple> output;

    for(auto tuple : table.getTuples())
    {
        if(query.evaluate(tuple, index))
        {
            output.push_back(tuple);
        }
    }

    //set up newvalues
    vector<string> newValues;
    for(auto fvPair : output.at(0).getFields())
    {
        newValues.push_back(fvPair.getField());
    }
    vector<string>::iterator it = newValues.begin();
    for(auto string : values)
    {
        it = find(newValues.begin(), newValues.end(), string);
        if(it != newValues.end())
        {
            newValues.erase(it);
        }
        
    }

    output = trim(output, newValues);

    return output;
}

vector<Tuple> nindexFunc(Relation table, OrNode query, set<int> indeces, vector<string> values)
{
    vector<Tuple> output;

    set<int>::iterator iter = indeces.begin();

    for(auto tuple : table.getTuples())
    {
        if(query.evaluate(tuple, (*iter)))
        {
            output.push_back(tuple);
        }
    }

    ++iter;
    for(; iter != indeces.end(); ++iter)
    {
        int i;
        for(i = 0; i < output.size(); i++)
        {
            if(!query.evaluate(output.at(i), (*iter)))
            {
                output.erase(output.begin() + i);
            }
        }
    }

    //set up newvalues
    vector<string> newValues;
    for(auto fvPair : output.at(0).getFields())
    {
        newValues.push_back(fvPair.getField());
    }
    vector<string>::iterator it = newValues.begin();
    for(auto string : values)
    {
        it = find(newValues.begin(), newValues.end(), string);
        if(it != newValues.end())
        {
            newValues.erase(it);
        }
        
    }

    output = trim(output, newValues);

    return output;
}

vector<Tuple> orderedindexFunc(Relation table, OrNode query, int index, vector<string> values)
{
    vector<Tuple> output;
    vector<Tuple> tablevalues = table.getTuples();
    tablevalues = sortFunc(tablevalues, index);

    for(auto tuple : tablevalues)
    {
        if(query.evaluate(tuple, index))
        {
            output.push_back(tuple);
        }
    }

    //set up newvalues
    vector<string> newValues;
    for(auto fvPair : output.at(0).getFields())
    {
        newValues.push_back(fvPair.getField());
    }
    vector<string>::iterator it = newValues.begin();
    for(auto string : values)
    {
        it = find(newValues.begin(), newValues.end(), string);
        if(it != newValues.end())
        {
            newValues.erase(it);
        }
        
    }

    output = trim(output, newValues);

    return output;
}

bool satisfiesJoin(Tuple tuple1, Tuple tuple2, int key1, int key2)
{
    if(tuple1.getValueAt(key1) == tuple2.getValueAt(key2))
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
        newTuple.addField(fvPair.getField(), fvPair.getValue());
    }
    for(FieldValuePair fvPair : two.getFields())
    {
        if(fvPair.getValue() != two.getValueAt(key))
        {
            newTuple.addField(fvPair.getField(), fvPair.getValue());
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
            tup1.getValueAt(key1), tup1));
    }
    for(auto tup2 : list2)
    {
        iter = hash.find(tup2.getValueAt(key2));
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
    list1 = sortFunc(list1, key1);
    list2 = sortFunc(list2, key2);

    vector<Tuple>::iterator iter1 = list1.begin();
    vector<Tuple>::iterator iter2 = list2.begin();

    while(iter1 != list1.end() && iter2 != list2.end())
    {
        if(satisfiesJoin((*iter1), (*iter2), key1, key2))
        {
            Tuple toAdd = joinTuples((*iter1), (*iter2), key2);
            output.push_back(toAdd);
            ++iter2;
        }
        else if((*iter1).compareTo((*iter2), key1, key2))
        {
        	if((*iter1).getFields().at(key1).getValue() == (*iter2).getFields().at(key2).getValue())
        	{
        		++iter2;
        	}
        	else
        	{
        		++iter1;
        	}
            
        }
        else
        {
            ++iter2;
        }
    }

    return output;
}

vector<Tuple> leftOuterJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2)
{
    vector<Tuple> output;

    for(auto tup1 : list1)
    {
        bool hadPair = false;
        for(auto tup2 : list2)
        {
            if(satisfiesJoin(tup1, tup2, key1, key2))
            {
                hadPair = true;
                Tuple toAdd = joinTuples(tup1, tup2, key2);
                output.push_back(toAdd);
            }
        }
        if(hadPair == false)
        {
            Tuple temp = list2.at(0);
            Tuple blank;
            for(auto fvPair : temp.getFields())
            {
                if(fvPair.getValue() == temp.getValueAt(key2))
                {
                    blank.addField(fvPair.getField(), fvPair.getValue());
                }
                else
                {
                    blank.addField("", "");
                }
            }
            Tuple toAdd = joinTuples(tup1, blank, key2);
            output.push_back(toAdd);
        }
        hadPair = false;
    }

    return output;
}

vector<Tuple> rightOuterJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2)
{
    vector<Tuple> output;

    for(auto tup1 : list2)
    {
        bool hadPair = false;
        for(auto tup2 : list1)
        {
            if(satisfiesJoin(tup1, tup2, key1, key2))
            {
                hadPair = true;
                Tuple toAdd = joinTuples(tup1, tup2, key2);
                output.push_back(toAdd);
            }
        }
        if(hadPair == false)
        {
            Tuple temp = list1.at(0);
            Tuple blank;
            for(auto fvPair : temp.getFields())
            {
                if(fvPair.getValue() == temp.getValueAt(key1))
                {
                    blank.addField(fvPair.getField(), fvPair.getValue());
                }
                else
                {
                    blank.addField("", "");
                }
            }
            Tuple toAdd = joinTuples(tup1, blank, key1);
            output.push_back(toAdd);
        }
        hadPair = false;
    }

    return output;
}




vector<Tuple> mergeFunc(vector<Tuple> left, vector<Tuple> right, int key1, int key2)
{
    vector<Tuple> output;

    while(left.size() != 0 && right.size() != 0)
    {
        if((*left.begin()).compareTo((*right.begin()), key1, key2))
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

vector<Tuple> sortFunc(vector<Tuple> list, int key)
{
    if(list.size() <= 1)
    {
        return list;
    }

    int mid = list.size() / 2;
    vector<Tuple>::iterator iter = list.begin() + mid;
    vector<Tuple> left(list.begin(), iter);
    vector<Tuple> right(iter, list.end());

    left = sortFunc(left, key);
    right = sortFunc(right, key);


    return mergeFunc(left, right, key, key);

}

vector<Tuple> unionFunc(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2)
{
    list1 = sortFunc(list1, key1);
    list2 = sortFunc(list2, key2);

    vector<Tuple> output;

    while(list1.size() != 0 && list2.size() != 0)
    {
        if((*list1.begin()).equals((*list2.begin())))
        {
            output.push_back((*list1.begin()));
            list1.erase(list1.begin());
            list2.erase(list2.begin());
        }
        else if((*list1.begin()).compareTo((*list2.begin()), key1, key2))
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

vector<Tuple> crossProduct(vector<Tuple> list1, vector<Tuple> list2)
{
    vector<Tuple> output;

    for(auto tup1 : list1)
    {
        for(auto tup2 : list2)
        {
            Tuple toAdd;
            for(auto fvPair : tup1.getFields())
            {
                toAdd.addField(fvPair.getField(), fvPair.getValue());
            }
            for(auto fvPair : tup2.getFields())
            {
                toAdd.addField(fvPair.getField(), fvPair.getValue());
            }
            output.push_back(toAdd);
        }
    }

    return output;
}

vector<Tuple> projection(vector<Tuple> list, vector<string> values)
{
    vector<Tuple> output;
    bool success = false;

    for(auto tuple : list)
    {
        for(auto string : values)
        {
            success = false;
            vector<FieldValuePair> search = tuple.getFields();
            for(auto fvPair : search)
            {
                if(string == fvPair.getField())
                {
                    success = true;
                    break;
                }
            }
            if(success != true)
            {
                break;
            }
        }
        
        if(success)
        {
            output.push_back(tuple);
        }
    }

    return output;
}

vector<Tuple> trim(vector<Tuple> list, vector<string> values)
{
    vector<Tuple> output;
    for(auto tuple : list)
    {
        for(auto string : values)
        {
            tuple.removeField(string);
        }
        output.push_back(tuple);
    }

    return output;
}

vector<Tuple> fullOuterJoin(vector<Tuple> list1, vector<Tuple> list2, int key1, int key2)
{
    vector<Tuple> output;

    vector<Tuple> left = leftOuterJoin(list1, list2, key1, key2);
    vector<Tuple> right = rightOuterJoin(list1, list2, key1, key2);
    left = sortFunc(left, key1);
    right = sortFunc(right, key2);

    output = unionFunc(left, right, key1, key2);

    return output;
}

void printTuples(vector<Tuple> list)
{
    for(auto tuple : list)
    {
        tuple.printTuple();

    }
}

#endif