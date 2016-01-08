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

void BuildUserTables()
{
	
}

RealPSet* UserFunction()
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

  Relation *orders = getRelationByName("orders");
  Relation *odetails = getRelationByName("odetails");

  InputNode *inA = new InputNode("orders", "ono", AFields, orders->getName(), "ono > 1000");
  InputNode *inB = new InputNode("odetails", "ono", BFields, odetails->getName(), "ono > 1000");
  inA->SetRelation(orders);
  inB->SetRelation(odetails);

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


/*Example Functions*/
void BuildExampleTables()
{
	//For zipcodes
 /* Relation *zipcodes = new Relation("zipcodes");
  zipcodes->addAttribute("zip", "number", true);
  zipcodes->addAttribute("city", "string", false);

  Tuple zip1;
  zip1.addField("zip", "67226");
  zip1.addField("city", "Wichita");
  zipcodes->addTuple(zip1);
  Tuple zip2;
  zip2.addField("zip", "60606");
  zip2.addField("city", "Fort Dodge");
  zipcodes->addTuple(zip2);
  Tuple zip3;
  zip3.addField("zip", "50302");
  zip3.addField("city", "Kansas City");
  zipcodes->addTuple(zip3);
  Tuple zip4;
  zip4.addField("zip", "54444");
  zip4.addField("city", "Columbia");
  zipcodes->addTuple(zip4);
  Tuple zip5;
  zip5.addField("zip", "66002");
  zip5.addField("city", "Liberal");
  zipcodes->addTuple(zip5);
  Tuple zip6;
  zip6.addField("zip", "61111");
  zip6.addField("city", "Fort Hays");
  zipcodes->addTuple(zip6);

  userRelations.push_back(zipcodes);

  //For employees
  Relation *employees = new Relation("employees");
  employees->addAttribute("eno", "number", true);
  employees->addAttribute("ename", "string", false);
  employees->addAttribute("zip", "number", true);
  employees->addAttribute("hdate", "string", false);

  Tuple emp1;
  emp1.addField("eno", "1000");
  emp1.addField("ename", "Jones");
  emp1.addField("zip", "67226");
  emp1.addField("hdate", "12-DEC-95");
  employees->addTuple(emp1);

  Tuple emp2;
  emp2.addField("eno", "1002");
  emp2.addField("ename", "Smith");
  emp2.addField("zip", "60606");
  emp2.addField("hdate", "01-JAN-92");
  employees->addTuple(emp2);

  Tuple emp3;
  emp3.addField("eno", "1002");
  emp3.addField("ename", "Brown");
  emp3.addField("zip", "50302");
  emp3.addField("hdate", "01-SEP-94");
  employees->addTuple(emp3);

  userRelations.push_back(employees);*/

  //For Parts
  Relation *parts = new Relation("parts");
  parts->addAttribute("pno", "number", true);
  parts->addAttribute("pname", "string", false);
  parts->addAttribute("qoh", "number", false);
  parts->addAttribute("price", "number", false);
  parts->addAttribute("olevel", "number", false);

  Tuple part1;
  part1.addField("pno", "10506");
  part1.addField("pname", "Land Before Time I");
  part1.addField("qoh", "200");
  part1.addField("price", "19.99");
  part1.addField("olevel", "20");
  parts->addTuple(part1);

  Tuple part2;
  part2.addField("pno", "10507");
  part2.addField("pname", "Land Before Time II");
  part2.addField("qoh", "156");
  part2.addField("price", "19.99");
  part2.addField("olevel", "20");
  parts->addTuple(part2);

  Tuple part3;
  part3.addField("pno", "10508");
  part3.addField("pname", "Land Before Time III");
  part3.addField("qoh", "190");
  part3.addField("price", "19.99");
  part3.addField("olevel", "20");
  parts->addTuple(part3);

  Tuple part4;
  part4.addField("pno", "10509");
  part4.addField("pname", "Land Before Time IV");
  part4.addField("qoh", "60");
  part4.addField("price", "19.99");
  part4.addField("olevel", "20");
  parts->addTuple(part4);

  Tuple part5;
  part5.addField("pno", "10601");
  part5.addField("pname", "Sleeping Beauty");
  part5.addField("qoh", "300");
  part5.addField("price", "24.99");
  part5.addField("olevel", "20");
  parts->addTuple(part5);

  Tuple part6;
  part6.addField("pno", "10701");
  part6.addField("pname", "When Harry Met Sally");
  part6.addField("qoh", "120");
  part6.addField("price", "19.99");
  part6.addField("olevel", "30");
  parts->addTuple(part6);

  Tuple part7;
  part7.addField("pno", "10800");
  part7.addField("pname", "Dirty Harry");
  part7.addField("qoh", "140");
  part7.addField("price", "14.99");
  part7.addField("olevel", "30");
  parts->addTuple(part7);

  Tuple part8;
  part8.addField("pno", "10900");
  part8.addField("pname", "Dr. Zhivago");
  part8.addField("qoh", "100");
  part8.addField("price", "24.99");
  part8.addField("olevel", "30");
  parts->addTuple(part8);

  userRelations.push_back(parts);

  /*//for Customers
  Relation *customers = new Relation("cumstomers");
  customers->addAttribute("cno", "number", true);
  customers->addAttribute("cname", "string", false);
  customers->addAttribute("street", "string", false);
  customers->addAttribute("zip", "number", true);
  customers->addAttribute("phone", "string", false);

  Tuple cust1;
  cust1.addField("cno", "1111");
  cust1.addField("cname", "Charles");
  cust1.addField("street", "123 Main St.");
  cust1.addField("zip", "67226");
  cust1.addField("phone", "316-636-5555");
  customers->addTuple(cust1);

  Tuple cust2;
  cust2.addField("cno", "2222");
  cust2.addField("cname", "Bertram");
  cust2.addField("street", "237 Ash Avenue");
  cust2.addField("zip", "67226");
  cust2.addField("phone", "316-689-5555");
  customers->addTuple(cust2);

  Tuple cust3;
  cust3.addField("cno", "3333");
  cust3.addField("cname", "Barbara");
  cust3.addField("street", "111 Inwood St.");
  cust3.addField("zip", "60606");
  cust3.addField("phone", "316-111-1234");
  customers->addTuple(cust3);

  userRelations.push_back(customers);*/

  //for orders
  Relation *orders = new Relation("orders");
  orders->addAttribute("ono", "number", true);
  orders->addAttribute("cno", "number", true);
  orders->addAttribute("eno", "number", true);
  orders->addAttribute("shipped date", "string", false);
  orders->addAttribute("received date", "string", false);

  Tuple order1;
  order1.addField("ono", "1020");
  order1.addField("cno", "1111");
  order1.addField("eno", "1000");
  order1.addField("shipped date", "10-DEC-94");
  order1.addField("received date", "12-DEC-94");
  orders->addTuple(order1);

  Tuple order2;
  order2.addField("ono", "1021");
  order2.addField("cno", "1111");
  order2.addField("eno", "1000");
  order2.addField("shipped date", "12-JAN-95");
  order2.addField("received date", "15-JAN-95");
  orders->addTuple(order2);

  Tuple order3;
  order3.addField("ono", "1022");
  order3.addField("cno", "2222");
  order3.addField("eno", "1001");
  order3.addField("shipped date", "13-FEB-95");
  order3.addField("received date", "20-FEB-95");
  orders->addTuple(order3);

  Tuple order4;
  order4.addField("ono", "1023");
  order4.addField("cno", "3333");
  order4.addField("eno", "1000");
  order4.addField("shipped date", "20-JUN-97");
  order4.addField("received date", "");
  orders->addTuple(order4);

  userRelations.push_back(orders);

  //for odetails
  Relation *odetails = new Relation("odetails");
  odetails->addAttribute("ono", "number", true);
  odetails->addAttribute("pno", "number", true);
  odetails->addAttribute("qty", "number", false);

  Tuple otail1;
  otail1.addField("ono", "1020");
  otail1.addField("pno", "10506");
  otail1.addField("qty", "1");
  odetails->addTuple(otail1);

  Tuple otail2;
  otail2.addField("ono", "1020");
  otail2.addField("pno", "10507");
  otail2.addField("qty", "1");
  odetails->addTuple(otail2);

  Tuple otail3;
  otail3.addField("ono", "1020");
  otail3.addField("pno", "10508");
  otail3.addField("qty", "2");
  odetails->addTuple(otail3);

  Tuple otail4;
  otail4.addField("ono", "1020");
  otail4.addField("pno", "10509");
  otail4.addField("qty", "3");
  odetails->addTuple(otail4);

  Tuple otail5;
  otail5.addField("ono", "1021");
  otail5.addField("pno", "10601");
  otail5.addField("qty", "4");
  odetails->addTuple(otail5);

  Tuple otail6;
  otail6.addField("ono", "1022");
  otail6.addField("pno", "10601");
  otail6.addField("qty", "1");
  odetails->addTuple(otail6);

  Tuple otail7;
  otail7.addField("ono", "1022");
  otail7.addField("pno", "10701");
  otail7.addField("qty", "1");
  odetails->addTuple(otail7);

  Tuple otail8;
  otail8.addField("ono", "1023");
  otail8.addField("pno", "10800");
  otail8.addField("qty", "1");
  odetails->addTuple(otail8);

  Tuple otail9;
  otail9.addField("ono", "1023");
  otail9.addField("pno", "10900");
  otail9.addField("qty", "1");
  odetails->addTuple(otail9);

  userRelations.push_back(odetails);
}

RealPSet* ExampleFunc()
{
  set<string> AFields;
  AFields.insert("ono");
  AFields.insert("cno");
  AFields.insert("eno");


  set<string> BFields;
  BFields.insert("ono");
  BFields.insert("pno");
  BFields.insert("qty");

  Relation *orders = getRelationByName("orders");
  Relation *odetails = getRelationByName("odetails");

  InputNode *inA = new InputNode("orders", "ono", AFields, orders->getName(), "ono > 1000");
  InputNode *inB = new InputNode("odetails", "ono", BFields, odetails->getName(), "ono > 1000");
  inA->SetRelation(orders);
  inB->SetRelation(odetails);

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

vector<Relation*> getUserRelations()
{
	return userRelations;
}

Relation* getRelationByName(string name)
{

  for(auto rel : getUserRelations())
  {
    if(name == rel->getName())
    {
      return rel;
    }
  }
  return getUserRelations().at(0);
}


#endif
