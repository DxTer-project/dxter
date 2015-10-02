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


#include "rqoTuple.h"

#if DORQO

bool Tuple::compareTo(Tuple comparator, int key)
{
    return (fields.at(key).getValue() < comparator.getFields().at(key).getValue());
}

bool Tuple::equals(Tuple comparator)
{
    vector<FieldValuePair>::iterator iter0 = fields.begin();
    vector<FieldValuePair>::iterator iter1 = comparator.getFields().begin();  
    for(; iter0 != fields.end() && iter1 != comparator.getFields().end(); ++iter0, ++iter1) 
    {
        if((*iter0).getField() != (*iter1).getField())
        {
            return false;
        }
        if((*iter0).getValue() != (*iter1).getValue())
        {
            return false;
        }
    }
    if(iter0 != fields.end() || iter1 != comparator.getFields().end())
    {
        return false;
    }
    return true;
}

#endif