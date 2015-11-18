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

#include "userInput.h"


#if DORQO

RealPSet* UserInput::ExampleFunc()
{
  set<string> AFields;
  AFields.insert("ono");
  AFields.insert("cno");
  AFields.insert("eno");
  AFields.insert("received");
  AFields.insert("shipped");

  set<string> BFields;
  BFields.insert("ono");
  BFields.insert("pno");
  BFields.insert("qty");


  Scan *inA = new Scan("orders", "ono", AFields, orders.getName(), "ono > 1000");
  Scan *inB = new Scan("odetails", "ono", BFields, odetails.getName(), "ono > 1000");
  inA->SetRelation(&orders);
  inB->SetRelation(&odetails);

  vector<string> joinFields0;
  joinFields0.push_back("ono");

  vector<string> joinFields1;
  joinFields1.push_back("ono");

  Join *join = new Join("ono", joinFields0, joinFields1);

  join->AddInput(inA, 0);
  join->AddInput(inB, 0);


  Poss *poss = new Poss(1, join);
  RealPSet *pset = new RealPSet(poss);

  return pset;
}



#endif
