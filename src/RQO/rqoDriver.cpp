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



#include "base.h"
#include "costs.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "transform.h"
#include "logging.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "rqoHelperNodes.h"
#include "rqoJoin.h"
#include "rqoProj.h"
#include "rqoSort.h"
#include "sortable.h"
#include "hJoin.h"
#include "rqoNode.h"
#include "functions.h"
#include "rqoRelation.h"
#include "rqoAttribute.h"
#include "rqoScan.h"
#include <sstream>
#include <unordered_map>
#include "userInput.h"

#if DORQO

//Database tables
Relation zipcodes("zipcodes");
Relation employees("employees");
Relation parts("parts");
Relation customers("customers");
Relation orders("orders");
Relation odetails("odetails");

unordered_map<string, vector<Tuple>> tuples;

UserInput user;

RealPSet* Example1();
RealPSet* Example2();
RealPSet* Example3();
RealPSet* Example4();
RealPSet* Example5();
void Example1Run();
void buildDatabases();
void parseCode(int bestAlg);
Relation* getRelation(string name);
FieldValue createQuery(string query);
int getKey(vector<Tuple> list, string fieldName);

typedef std::chrono::time_point<std::chrono::system_clock> AccurateTime;

double difftime(AccurateTime &end, AccurateTime &start)
{
  return (std::chrono::duration_cast<std::chrono::milliseconds>(end-start)).count()/1000.0;
}


void AddTrans()
{
  //  Universe::AddTrans(Projection::GetClass(), new RemoveExtraProjection, RQOPHASE);
 // Universe::AddTrans(Join::GetClass(), new SwapNodes(0,Join::GetClass()), RQOPHASE);
 // Universe::AddTrans(HJoin::GetClass(), new SwapNodes(0,HJoin::GetClass()), RQOPHASE);
  //Universe::AddTrans(Join::GetClass(), new SwapNodes(1,Join::GetClass()), RQOPHASE);
  //Universe::AddTrans(HJoin::GetClass(), new SwapNodes(1,HJoin::GetClass()), RQOPHASE);
  Universe::AddTrans(Join::GetClass(), new JoinToHash, RQOPHASE);
  Universe::AddTrans(Join::GetClass(), new JoinToNested, RQOPHASE);
  Universe::AddTrans(Join::GetClass(), new JoinToMerge, RQOPHASE);
}

void AddSimplifiers()
{ 
  //Universe::AddTrans(Projection::GetClass(), new RemoveExtraProjection, SIMP);
  
}

void Usage()
{
  cout << "./driver \n";
  //  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> example1\n";
  cout <<"         2  -> example2\n";
  cout <<"         3  -> example3\n";
  cout <<"         4  -> example4\n";
  cout <<"         5  -> example5\n";
}

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_nested(true);
  omp_init_lock(&RealPSet::m_lock);
#endif
  LOG_START("tensors");
  //  PrintType printType = CODE;
  RealPSet* (*algFunc)();
  //  GraphNum whichGraph = 0;
  int algNum;
  string fileName;
  buildDatabases();

  if(argc < 2) {
    Usage();
    return 0;
  }
  else {
    algNum = atoi(argv[1]);
    switch(algNum) {
    case(1):
      algFunc = Example1;
      break;
    case(2):
      algFunc = Example2;
      break;
    case(3):
      algFunc = Example3;
      break;
    case(4):
      algFunc = Example4;
      break;
    case(5):
      algFunc = Example5;
      break;
    /*case(6):
      algFunc = user.ExampleFunc;
      break;*/
    default:
      Usage();
      return 0;
    }
  }
  AddTrans();
  AddSimplifiers();

  Universe uni;
  AccurateTime start, start2, end;
  uni.PrintStats();

  if (algNum==0) {
    throw;
  }
  else {
    RealPSet *startSet = algFunc();
    uni.Init(startSet);
    uni.Prop();
    GraphIter graphIter(startSet->m_posses.begin()->second);
    Cost flopCost = graphIter.EvalAndSetBest();
    cout << "*****FLOPS = " << setprecision(15) << flopCost << endl;
    start = std::chrono::system_clock::now();
  }


#if DORQOPHASE
  if (CurrPhase == RQOPHASE) {
    start2 = std::chrono::system_clock::now();
    cout << "Expanding phase 1\n";
    uni.Expand(-1, RQOPHASE, NULL);
    end = std::chrono::system_clock::now();
    cout << "First phase took " << difftime(end,start2) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    start2 = std::chrono::system_clock::now();
    uni.Prop();
    end = std::chrono::system_clock::now();
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

  end = std::chrono::system_clock::now();
  cout << "Full expansion took " << difftime(end,start) << " seconds\n";
cout << "Left with " << uni.TotalCount() << " algorithms\n";
  cout.flush();
  
int best = 0;

#if 1
  best = uni.PrintAll(algNum);
#else
  uni.PrintBest();
#endif
  //cout << algNum << endl;
  //cout << best << endl;
  LOG_END();

  
 /* zipcodes.printTable();
  employees.printTable();
  parts.printTable();
  customers.printTable();
  orders.printTable();
  odetails.printTable();*/
  //Example1Run();
  parseCode(best);
  return 0;
}


void parseCode(int bestAlg)
{
  unordered_map<string, vector<Tuple>>::iterator iter;
  string line;
  ifstream infile;
  bool uniqueFound = false;
  infile.open ("codeOutput.txt");
  if(!infile)
  {
    cout << "Failed to find code file." << endl;
    infile.close();
    return;
  }
  stringstream test;
  //test << "\tUnique Num: " << bestAlg;

//get to algorithm we need
  while(!uniqueFound)
  {
    getline(infile, line);
    
    if(line == test.str())
    {
      break;
    }
    else if(infile.eof())
    {
      break;
    }
  }
  while(line != "//------------------------------------//")
  {
    getline(infile, line);
  }
  getline(infile, line);
  string last;
  while(line != "//------------------------------------//")
  {
    getline(infile, line);
    //begin the parsing
    if(line == "//------------------------------------//")
    {
      //done
      break;
    }
    else if(line != "")
    {
      //cout << line << endl;
      int i = 0;
      int j = 0;
      //get name of variable created
      i = line.find(" ");
      string returnName = line.substr(j, i);
      last = returnName;
      //cout << returnName << endl;
      i += 3;
      j = i;
      //get name of function we're using
      i = line.find("(");
      string funcName = line.substr(j, i - j);
      //cout << funcName << endl;
      j = i + 1;
      //if scan
      if(funcName == "scanFunc")
      {
        i = line.find(",");
        string relationName = line.substr(j, i - j);
        j = i + 1;
        i = line.find(")");
        string query = line.substr(j, i - j);
        Relation *tempRel = getRelation(relationName);
        FieldValue tempOr = createQuery(query);
        vector<Tuple> tempTup = scanFunc((*tempRel), tempOr);
        //save output
        tuples.insert(pair<string, vector<Tuple>>(
            returnName, tempTup));
      }
      else if(funcName == "indexFunc" || funcName == "orderedindexFunc")
      {
        i = line.find(",");
        string relationName = line.substr(j, i - j);
        j = i + 1;
        i = line.find(",", j);
        string query = line.substr(j, i - j);
        j = i + 1;
        i = line.find(")", j);
        string index = line.substr(j, i - j);
        Relation *tempRel = getRelation(relationName);
        FieldValue tempOr = createQuery(query);
        int key = getKey(tempRel->getTuples(), index);
        vector<Tuple> tempTup;
        if(funcName == "indexFunc")
        {
          tempTup = indexFunc((*tempRel), tempOr, key);
        }
        else
        {
          tempTup = orderedindexFunc((*tempRel), tempOr, key);
        }
        //save output
        tuples.insert(pair<string, vector<Tuple>>(
            returnName, tempTup));
      }
      else if(funcName == "nindexFunc")
      {
        i = line.find(",");
        string relationName = line.substr(j, i - j);
        j = i + 1;
        i = line.find(",", j);
        string query = line.substr(j, i - j);
        j = i + 2;

        Relation *tempRel = getRelation(relationName);
        FieldValue tempOr = createQuery(query);
        set<int> indeces;
        while(line.find("]") != j)
        {
          i = line.find(",", j);
          if(i == string::npos)
          {
            i = line.find("]");
          }
          string index = line.substr(j, i - j);
          int key = getKey(tempRel->getTuples(), index);
          indeces.insert(key);
          j = i + 1;
        }
        
        
        vector<Tuple> tempTup;
        tempTup = nindexFunc((*tempRel), tempOr, indeces);


        tuples.insert(pair<string, vector<Tuple>>(
            returnName, tempTup));
      }
      else if(funcName == "nestedJoin" 
        || funcName == "mergeJoin" 
        || funcName == "hashJoin")
      {
        i = line.find(",") + 1;
        j = i;
        i = line.find(",", j);
        string set1 = line.substr(j, i - j);
        j = i + 1;
        i = line.find(",", j);
        string set2 = line.substr(j, i - j);

        iter = tuples.find(set1);
        vector<Tuple> left = iter->second;
        iter = tuples.find(set2);
        vector<Tuple> right = iter->second;

        j = line.find(".", i) + 1;
        i = line.find(",", j);
        string field1 = line.substr(j, i - j);
        j = line.find(".", i) + 1;
        i = line.find(")");
        string field2 = line.substr(j, i - j);

        int key1 = getKey(left, field1);
        int key2 = getKey(right, field2);
        vector<Tuple> tempTup;

        if(funcName == "nestedJoin")
        {
          tempTup = nestedJoin(left, right, key1, key2);
        }
        else if(funcName == "mergeJoin")
        {
          tempTup = mergeJoin(left, right, key1, key2);
        }
        else
        {
          tempTup = hashJoin(left, right, key1, key2);
        }
        

        tuples.insert(pair<string, vector<Tuple>>(
            returnName, tempTup));

      }
      else if(funcName == "sortFunc")
      {
        i = line.find(",", j);
        string set = line.substr(j, i - j);
        j = i + 1;
        i = line.find(")", j);
        string sortBy = line.substr(j, i - j);

        iter = tuples.find(set);
        vector<Tuple> list = iter->second;
        int key = getKey(list, sortBy);

        vector<Tuple> tempTup = sortFunc(list, key);

        tuples.insert(pair<string, vector<Tuple>>(
            returnName, tempTup));
      }
      else if(funcName == "unionFunc")
      {
        i = line.find(",") + 1;
        j = i;
        i = line.find(",", j);
        string set1 = line.substr(j, i - j);
        j = i + 1;
        i = line.find(",", j);
        string set2 = line.substr(j, i - j);

        iter = tuples.find(set1);
        vector<Tuple> left = iter->second;
        iter = tuples.find(set2);
        vector<Tuple> right = iter->second;

        j = line.find(".", i) + 1;
        i = line.find(",", j);
        string field1 = line.substr(j, i - j);
        j = line.find(".", i) + 1;
        i = line.find(")");
        string field2 = line.substr(j, i - j);

        int key1 = getKey(left, field1);
        int key2 = getKey(right, field2);
        vector<Tuple> tempTup = unionFunc(left, right, key1, key2);

        tuples.insert(pair<string, vector<Tuple>>(
            returnName, tempTup));
      }
      else if(funcName == "crossProduct")
      {
        i = line.find(",", j);
        string set1 = line.substr(j, i - j);
        j = i + 1;
        i = line.find(",", j);
        string set2 = line.substr(j, i - j);

        iter = tuples.find(set1);
        vector<Tuple> left = iter->second;
        iter = tuples.find(set2);
        vector<Tuple> right = iter->second;

        vector<Tuple> tempTup = crossProduct(left, right);

        tuples.insert(pair<string, vector<Tuple>>(
            returnName, tempTup));
      }
      else if(funcName == "projection")
      {
        i = line.find(",") + 1;
        j = i;
        i = line.find(",", j);
        string set = line.substr(j, i - j);

        iter = tuples.find(set);
        vector<Tuple> left = iter->second;

        j = i + 1;
        vector<string> values;
        while(line.find(")") != j)
        {
          j = line.find(".", i) + 1;
          i = line.find(",", j);
          string tempVal = line.substr(j, i - j);
          values.push_back(tempVal);
          j = i + 1;
        }

        vector<Tuple> tempTup = projection(left, values);

        tuples.insert(pair<string, vector<Tuple>>(
            returnName, tempTup));
      }
      else
      {
        cout << "invalid function call" << endl;
        throw;
      }

    }

    
  }

  iter = tuples.find(last);
    if(iter != tuples.end())
    {
      vector<Tuple> final = iter->second;
      for(auto t : final)
      {
        t.printTuple();
      }
    }
    
  infile.close();
}

int getKey(vector<Tuple> list, string fieldName)
{
  vector<FieldValuePair> temp = list.at(1).fields;

  int i = 0;
  for(auto fvPair : temp)
  {
    if(fvPair.getField() == fieldName)
    {
      break;
    }
    i++;
  }
  return i;
}

Relation* getRelation(string name)
{
  Relation *match;


  if(name == zipcodes.getName())
  {
    match = &zipcodes;
  }
  else if(name == "employees")
  {
    match = &employees;
  }
  else if(name == "parts")
  {
    match = &parts;
  }
  else if(name == "customers")
  {
    match = &customers;
  }
  else if(name == "orders")
  {
    match = &orders;
  }
  else if(name == "odetails")
  {
    match = &odetails;
  }
  return match;
}

FieldValue createQuery(string query)
{
  string relation;
  int i = 0;
  int j = 0;
  i = query.find(" ");
  ++i;
  j = query.find(" ", i);
  relation = query.substr(i, j - i);
  ++j;
  string value = query.substr(j);
  FieldValue ret (relation, value);
  //OrNode blah;
 // AndNode yep;
  //yep.addClause(&ret);
  //blah.addAnd(&yep);
  return ret;
}

//void runGraphCode()

RealPSet* Example1()
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

/*void Example1Run()
{
  cout << "Beginning Functions\n\n";

  //create tuples
  vector<Tuple> inA;
  vector<Tuple> inB;
  vector<Tuple> inC;

//name, bank, SSN
    FieldValuePair aTemp11("x", "Jonathon");
    FieldValuePair aTemp12("y", "Wells Fargo");
    FieldValuePair aTemp13("z", "XXX-XX-XXXX");
    Tuple atemp1;
    atemp1.addField(aTemp11);
    atemp1.addField(aTemp12);
    atemp1.addField(aTemp13);
    inA.push_back(atemp1);
    FieldValuePair aTemp21("x", "Mary");
    FieldValuePair aTemp22("y", "Wells Fargo");
    FieldValuePair aTemp23("z", "XXX-XX-XXXX");
    Tuple atemp2;
    atemp2.addField(aTemp21);
    atemp2.addField(aTemp22);
    atemp2.addField(aTemp23);
    inA.push_back(atemp2);
    FieldValuePair aTemp31("x", "Alfred");
    FieldValuePair aTemp32("y", "Bank of America");
    FieldValuePair aTemp33("z", "XXX-XX-XXXX");
    Tuple atemp3;
    atemp3.addField(aTemp31);
    atemp3.addField(aTemp32);
    atemp3.addField(aTemp33);
    inA.push_back(atemp3);

//name, age, classification
    FieldValuePair bTemp11("u", "Jonathon");
    FieldValuePair bTemp12("v", "22");
    FieldValuePair bTemp13("w", "Junior Member");
    Tuple btemp1;
    btemp1.addField(bTemp11);
    btemp1.addField(bTemp12);
    btemp1.addField(bTemp13);
    inB.push_back(btemp1);
    FieldValuePair bTemp21("u", "Alfred");
    FieldValuePair bTemp22("v", "65");
    FieldValuePair bTemp23("w", "Senior Member");
    Tuple btemp2;
    btemp2.addField(bTemp21);
    btemp2.addField(bTemp22);
    btemp2.addField(bTemp23);
    inB.push_back(btemp2);
    FieldValuePair bTemp31("u", "Jorge");
    FieldValuePair bTemp32("v", "18");
    FieldValuePair bTemp33("w", "Junior Member");
    Tuple btemp3;
    btemp3.addField(bTemp31);
    btemp3.addField(bTemp32);
    btemp3.addField(bTemp33);
    inB.push_back(btemp3);

//state, city, bank
    FieldValuePair cTemp11("a", "Colorado");
    FieldValuePair cTemp12("b", "Denver");
    FieldValuePair cTemp13("c", "Wells Fargo");
    Tuple ctemp1;
    ctemp1.addField(cTemp11);
    ctemp1.addField(cTemp12);
    ctemp1.addField(cTemp13);
    inC.push_back(ctemp1);
    FieldValuePair cTemp21("a", "Tennessee");
    FieldValuePair cTemp22("b", "Memphis");
    FieldValuePair cTemp23("c", "Bank of America");
    Tuple ctemp2;
    ctemp2.addField(cTemp21);
    ctemp2.addField(cTemp22);
    ctemp2.addField(cTemp23);
    inC.push_back(ctemp2);
   

  //call functions
  cout << "Input : " << endl;
  printTuples(inA);
  printTuples(inB);
  printTuples(inC);

  vector<Tuple> fromJoin1 = hashJoin(inA, inB, 0, 0);

  cout << "hJoin1 output : " << endl;
  printTuples(fromJoin1);
  vector<Tuple> fromJoin2 = hashJoin(fromJoin1, inC, 1, 2);
  cout << "hjoin2 output : " << endl;
  printTuples(fromJoin2);

  vector<string> projValues;
  projValues.push_back("b");
  projValues.push_back("y");
  projValues.push_back("x");

  vector<Tuple> output = projection(fromJoin2, projValues);

  //show results
  cout << "after projection : " << endl;
  printTuples(fromJoin2);
}*/

RealPSet* Example2()
{
  set<string> AFields;
  AFields.insert("x");
  AFields.insert("y");
  AFields.insert("z");

  set<string> CFields;
  CFields.insert("a");
  CFields.insert("b");
  CFields.insert("c");

  InputNode *inA = new InputNode("A", "x", AFields, "fileName", "query");
  InputNode *inC = new InputNode("C", "a", CFields, "fileName", "query");

  set<string> projFields1;
  projFields1.insert("x");  
  projFields1.insert("y");
 

  set<string> projFields2;
  projFields2.insert("b");
  projFields2.insert("a");

  Projection *proj1 = new Projection("x", projFields1);

  Projection *proj2 = new Projection("b", projFields2);

  proj1->AddInput(inA, 0);
  proj2->AddInput(inC, 0);

  vector<string> joinFields0;
  joinFields0.push_back("x");
  joinFields0.push_back("y");

  vector<string> joinFields1;
  joinFields1.push_back("a");
  joinFields1.push_back("b");

  Join *join = new Join("x", joinFields0, joinFields1);

  join->AddInput(proj1, 0);
  join->AddInput(proj2, 0);

  Poss *poss = new Poss(1, join);
  RealPSet *pset = new RealPSet(poss);
  return pset;
}

RealPSet* Example3()
{
  set<string> AFields;
  AFields.insert("x");
  AFields.insert("y");
  AFields.insert("z");

  InputNode *inA = new InputNode("A", "x", AFields, "fileName", "query");

  set<string> projFields1;
  projFields1.insert("x");  
  projFields1.insert("y");
 

  set<string> projFields2;
  projFields2.insert("x");


  Projection *proj1 = new Projection("x", projFields1);

  Projection *proj2 = new Projection("x", projFields2);

  proj1->AddInput(inA, 0);

  proj2->AddInput(proj1, 0);

  Poss *poss = new Poss(1, proj2);
  RealPSet *pset = new RealPSet(poss);
  return pset;
}

RealPSet* Example4()
{
  set<string> AFields;
  AFields.insert("x");
  AFields.insert("y");
  AFields.insert("z");

  set<string> BFields;
  BFields.insert("u");
  BFields.insert("v");
  BFields.insert("w");

  set<string> CFields;
  CFields.insert("a");
  CFields.insert("b");
  CFields.insert("c");

  InputNode *inA = new InputNode("A", "x", AFields, "fileName", "query");
  InputNode *inB = new InputNode("B", "u", BFields, "fileName", "query");
  InputNode *inC = new InputNode("C", "a", CFields, "fileName", "query");

  vector<string> joinFields0;
  joinFields0.push_back("x");

  vector<string> joinFields1;
  joinFields1.push_back("u");

  Join *join = new Join("u", joinFields0, joinFields1);

  join->AddInput(inA, 0);
  join->AddInput(inB, 0);


  vector<string> joinFields2;
  joinFields2.push_back("x");

  vector<string> joinFields3;
  joinFields3.push_back("a");

  Join *join2 = new Join("x", joinFields2, joinFields3);

  join2->AddInput(join, 0);
  join2->AddInput(inC, 0);


  Poss *poss = new Poss(1, join2);
  RealPSet *pset = new RealPSet(poss);
  return pset;
}

RealPSet* Example5()
{
  set<string> AFields;
  AFields.insert("x");
  AFields.insert("y");
  AFields.insert("z");

  set<string> BFields;
  BFields.insert("u");
  BFields.insert("v");
  BFields.insert("w");

  set<string> CFields;
  CFields.insert("a");
  CFields.insert("b");
  CFields.insert("c");

  InputNode *inA = new InputNode("A", "x", AFields, "fileName", "query");
  InputNode *inB = new InputNode("B", "u", BFields, "fileName", "query");
  InputNode *inC = new InputNode("C", "a", CFields, "fileName", "query");

  vector<string> joinFields0;
  joinFields0.push_back("x");

  vector<string> joinFields1;
  joinFields1.push_back("u");

  Join *join = new Join("u", joinFields0, joinFields1);

  join->AddInput(inA, 0);
  join->AddInput(inB, 0);


  vector<string> joinFields2;
  joinFields2.push_back("y");

  vector<string> joinFields3;
  joinFields3.push_back("a");

  Join *join2 = new Join("a", joinFields2, joinFields3);

  join2->AddInput(join, 0);
  join2->AddInput(inC, 0);


  Poss *poss = new Poss(1, join2);
  RealPSet *pset = new RealPSet(poss);
  return pset;
}

//This method will create and populate relation tables for the program to use
void buildDatabases()
{
  //For zipcodes
  zipcodes.addAttribute("zip", "number", true);
  zipcodes.addAttribute("city", "string", false);

  Tuple zip1;
  zip1.addField("zip", "67226");
  zip1.addField("city", "Wichita");
  zipcodes.addTuple(zip1);
  Tuple zip2;
  zip2.addField("zip", "60606");
  zip2.addField("city", "Fort Dodge");
  zipcodes.addTuple(zip2);
  Tuple zip3;
  zip3.addField("zip", "50302");
  zip3.addField("city", "Kansas City");
  zipcodes.addTuple(zip3);
  Tuple zip4;
  zip4.addField("zip", "54444");
  zip4.addField("city", "Columbia");
  zipcodes.addTuple(zip4);
  Tuple zip5;
  zip5.addField("zip", "66002");
  zip5.addField("city", "Liberal");
  zipcodes.addTuple(zip5);
  Tuple zip6;
  zip6.addField("zip", "61111");
  zip6.addField("city", "Fort Hays");
  zipcodes.addTuple(zip6);

  //For employees
  employees.addAttribute("eno", "number", true);
  employees.addAttribute("ename", "string", false);
  employees.addAttribute("zip", "number", true);
  employees.addAttribute("hdate", "string", false);

  Tuple emp1;
  emp1.addField("eno", "1000");
  emp1.addField("ename", "Jones");
  emp1.addField("zip", "67226");
  emp1.addField("hdate", "12-DEC-95");
  employees.addTuple(emp1);

  Tuple emp2;
  emp2.addField("eno", "1002");
  emp2.addField("ename", "Smith");
  emp2.addField("zip", "60606");
  emp2.addField("hdate", "01-JAN-92");
  employees.addTuple(emp2);

  Tuple emp3;
  emp3.addField("eno", "1002");
  emp3.addField("ename", "Brown");
  emp3.addField("zip", "50302");
  emp3.addField("hdate", "01-SEP-94");
  employees.addTuple(emp3);

  //For Parts
  parts.addAttribute("pno", "number", true);
  parts.addAttribute("pname", "string", false);
  parts.addAttribute("qoh", "number", false);
  parts.addAttribute("price", "number", false);
  parts.addAttribute("olevel", "number", false);

  Tuple part1;
  part1.addField("pno", "10506");
  part1.addField("pname", "Land Before Time I");
  part1.addField("qoh", "200");
  part1.addField("price", "19.99");
  part1.addField("olevel", "20");
  parts.addTuple(part1);

  Tuple part2;
  part2.addField("pno", "10507");
  part2.addField("pname", "Land Before Time II");
  part2.addField("qoh", "156");
  part2.addField("price", "19.99");
  part2.addField("olevel", "20");
  parts.addTuple(part2);

  Tuple part3;
  part3.addField("pno", "10508");
  part3.addField("pname", "Land Before Time III");
  part3.addField("qoh", "190");
  part3.addField("price", "19.99");
  part3.addField("olevel", "20");
  parts.addTuple(part3);

  Tuple part4;
  part4.addField("pno", "10509");
  part4.addField("pname", "Land Before Time IV");
  part4.addField("qoh", "60");
  part4.addField("price", "19.99");
  part4.addField("olevel", "20");
  parts.addTuple(part4);

  Tuple part5;
  part5.addField("pno", "10601");
  part5.addField("pname", "Sleeping Beauty");
  part5.addField("qoh", "300");
  part5.addField("price", "24.99");
  part5.addField("olevel", "20");
  parts.addTuple(part5);

  Tuple part6;
  part6.addField("pno", "10701");
  part6.addField("pname", "When Harry Met Sally");
  part6.addField("qoh", "120");
  part6.addField("price", "19.99");
  part6.addField("olevel", "30");
  parts.addTuple(part6);

  Tuple part7;
  part7.addField("pno", "10800");
  part7.addField("pname", "Dirty Harry");
  part7.addField("qoh", "140");
  part7.addField("price", "14.99");
  part7.addField("olevel", "30");
  parts.addTuple(part7);

  Tuple part8;
  part8.addField("pno", "10900");
  part8.addField("pname", "Dr. Zhivago");
  part8.addField("qoh", "100");
  part8.addField("price", "24.99");
  part8.addField("olevel", "30");
  parts.addTuple(part8);


  //for Customers
  customers.addAttribute("cno", "number", true);
  customers.addAttribute("cname", "string", false);
  customers.addAttribute("street", "string", false);
  customers.addAttribute("zip", "number", true);
  customers.addAttribute("phone", "string", false);

  Tuple cust1;
  cust1.addField("cno", "1111");
  cust1.addField("cname", "Charles");
  cust1.addField("street", "123 Main St.");
  cust1.addField("zip", "67226");
  cust1.addField("phone", "316-636-5555");
  customers.addTuple(cust1);

  Tuple cust2;
  cust2.addField("cno", "2222");
  cust2.addField("cname", "Bertram");
  cust2.addField("street", "237 Ash Avenue");
  cust2.addField("zip", "67226");
  cust2.addField("phone", "316-689-5555");
  customers.addTuple(cust2);

  Tuple cust3;
  cust3.addField("cno", "3333");
  cust3.addField("cname", "Barbara");
  cust3.addField("street", "111 Inwood St.");
  cust3.addField("zip", "60606");
  cust3.addField("phone", "316-111-1234");
  customers.addTuple(cust3);

  //for orders
  orders.addAttribute("ono", "number", true);
  orders.addAttribute("cno", "number", true);
  orders.addAttribute("eno", "number", true);
  orders.addAttribute("shipped date", "string", false);
  orders.addAttribute("received date", "string", false);

  Tuple order1;
  order1.addField("ono", "1020");
  order1.addField("cno", "1111");
  order1.addField("eno", "1000");
  order1.addField("shipped date", "10-DEC-94");
  order1.addField("received date", "12-DEC-94");
  orders.addTuple(order1);

  Tuple order2;
  order2.addField("ono", "1021");
  order2.addField("cno", "1111");
  order2.addField("eno", "1000");
  order2.addField("shipped date", "12-JAN-95");
  order2.addField("received date", "15-JAN-95");
  orders.addTuple(order2);

  Tuple order3;
  order3.addField("ono", "1022");
  order3.addField("cno", "2222");
  order3.addField("eno", "1001");
  order3.addField("shipped date", "13-FEB-95");
  order3.addField("received date", "20-FEB-95");
  orders.addTuple(order3);

  Tuple order4;
  order4.addField("ono", "1023");
  order4.addField("cno", "3333");
  order4.addField("eno", "1000");
  order4.addField("shipped date", "20-JUN-97");
  order4.addField("received date", "");
  orders.addTuple(order4);

  //for odetails
  odetails.addAttribute("ono", "number", true);
  odetails.addAttribute("pno", "number", true);
  odetails.addAttribute("qty", "number", false);

  Tuple otail1;
  otail1.addField("ono", "1020");
  otail1.addField("pno", "10506");
  otail1.addField("qty", "1");
  odetails.addTuple(otail1);

  Tuple otail2;
  otail2.addField("ono", "1020");
  otail2.addField("pno", "10507");
  otail2.addField("qty", "1");
  odetails.addTuple(otail2);

  Tuple otail3;
  otail3.addField("ono", "1020");
  otail3.addField("pno", "10508");
  otail3.addField("qty", "2");
  odetails.addTuple(otail3);

  Tuple otail4;
  otail4.addField("ono", "1020");
  otail4.addField("pno", "10509");
  otail4.addField("qty", "3");
  odetails.addTuple(otail4);

  Tuple otail5;
  otail5.addField("ono", "1021");
  otail5.addField("pno", "10601");
  otail5.addField("qty", "4");
  odetails.addTuple(otail5);

  Tuple otail6;
  otail6.addField("ono", "1022");
  otail6.addField("pno", "10601");
  otail6.addField("qty", "1");
  odetails.addTuple(otail6);

  Tuple otail7;
  otail7.addField("ono", "1022");
  otail7.addField("pno", "10701");
  otail7.addField("qty", "1");
  odetails.addTuple(otail7);

  Tuple otail8;
  otail8.addField("ono", "1023");
  otail8.addField("pno", "10800");
  otail8.addField("qty", "1");
  odetails.addTuple(otail8);

  Tuple otail9;
  otail9.addField("ono", "1023");
  otail9.addField("pno", "10900");
  otail9.addField("qty", "1");
  odetails.addTuple(otail9);

}


#endif


