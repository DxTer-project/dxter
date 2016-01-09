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
#include "userInput.h"
#include "functions.h"
#include "rqoRelation.h"
#include "rqoAttribute.h"
#include "rqoScan.h"
#include "rqoTrim.h"
#include <sstream>
#include <unordered_map>


#if DORQO

vector<Relation*> userRelations;

unordered_map<string, vector<Tuple>> tuples;

RealPSet* UserFunc();
RealPSet* Example();
void parseCode(int bestAlg);
int getRelationSpot(string name);
Relation* getRelation(string name);
OrNode* createQuery(string query);
AndNode* findAnd(string chunk);
FieldValue* findClause(string query);
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
  Universe::AddTrans(InputNode::GetClass(), new InputToScan, RQOPHASE);
  Universe::AddTrans(InputNode::GetClass(), new InputToIndex, RQOPHASE);
  Universe::AddTrans(InputNode::GetClass(), new InputToNIndex, RQOPHASE);
  Universe::AddTrans(InputNode::GetClass(), new InputToOrdered, RQOPHASE);
}

void AddSimplifiers()
{ 
  //Universe::AddTrans(Projection::GetClass(), new RemoveExtraProjection, SIMP);
  
}

void Usage()
{
  cout << "./driver \n";
  //  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Your User Function\n";
  cout <<"         2  -> The Example Function\n";
}

int main(int argc, const char* argv[])
{
  try
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

  if(argc < 2) {
    Usage();
    return 0;
  }
  else {
    algNum = atoi(argv[1]);
    switch(algNum) {
    case(1):
      algFunc = UserFunc;
      BuildUserTables();
      break;
    case(2):
      algFunc = Example;
      BuildExampleTables();
      break;
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
  LOG_END();

  parseCode(best);
}
catch(...)
{
  cout << "An Error occured. Please check to make sure you have not made any typos"
    << " and that all of your code statements follow proper syntax." << endl;
}
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
        i = line.find(",", j);
        string sortBy = line.substr(j, i - j);
        j = i + 1;
        i = line.find(",", j);
        string query = line.substr(j, i - j);

        j = i + 1;
        vector<string> values;
        while(line.find(")") > j)
        {
          i = line.find(",", j);
          if(i == -1)
          {
            i = line.find(")", j);
          }
          string tempVal = line.substr(j, i - j);
          values.push_back(tempVal);
          j = i + 1;
        }


        Relation *tempRel = getRelation(relationName);
        OrNode *tempOr = createQuery(query);
        int sortOn = getKey(tempRel->getTuples(), sortBy);
        vector<Tuple> tempTup = scanFunc((*tempRel), tempOr, values);
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
        i = line.find(",", j);
        string ind = line.substr(j, i - j);
        int index = atoi(ind.c_str());

        j = i + 1;
        vector<string> values;
        while(line.find(")") > j)
        {
          i = line.find(",", j);
          if(i == -1)
          {
            i = line.find(")", j);
          }
          string tempVal = line.substr(j, i - j);
          values.push_back(tempVal);
          j = i + 1;
        }

        Relation *tempRel = getRelation(relationName);
        OrNode *tempOr = createQuery(query);
        int key = index;
        vector<Tuple> tempTup;
        if(funcName == "indexFunc")
        {
          tempTup = indexFunc((*tempRel), tempOr, key, values);
        }
        else
        {
          tempTup = orderedindexFunc((*tempRel), tempOr, key, values);
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
        OrNode *tempOr = createQuery(query);
        set<int> indeces;
        while(line.find("]") != j)
        {
          i = line.find(",", j);
          if(i > line.find("]"))
          {
            i = line.find("]");
          }
          string index = line.substr(j, i - j);
          int key = atoi(index.c_str());
          indeces.insert(key);
          j = i + 1;
        }
        
        vector<string> values;
        while(line.find(")") > j)
        {
          i = line.find(",", j);
          if(i == -1)
          {
            i = line.find(")", j);
          }
          string tempVal = line.substr(j, i - j);
          values.push_back(tempVal);
          j = i + 1;
        }
        
        vector<Tuple> tempTup;
        tempTup = nindexFunc((*tempRel), tempOr, indeces, values);


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
      else if(funcName == "projection" || funcName == "trim")
      {
        i = line.find(",") + 1;
        j = i;
        i = line.find(",", j);
        string set = line.substr(j, i - j);

        iter = tuples.find(set);
        vector<Tuple> left = iter->second;

        j = i + 1;
        vector<string> values;
        while(line.find(")") > j)
        {
          j = line.find(".", i) + 1;
          i = line.find(",", j);
          if(i == -1)
          {
            i = line.find(")", j);
          }
          string tempVal = line.substr(j, i - j);
          values.push_back(tempVal);
          j = i + 1;
        }

        vector<Tuple> tempTup;
        if(funcName == "projection")
        {
          tempTup = projection(left, values);
        }
        else
        {
          tempTup = trim(left, values);
        }
        

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

int getRelationSpot(string name)
{
  int spot;

  for(int i = 0; i < userRelations.size(); i++)
  {
    if(name == userRelations.at(i)->getName())
    {
      spot = i;
      break;
    }
  }
  return spot;

}

Relation* getRelation(string name)
{

  return userRelations.at(getRelationSpot(name));
}

OrNode* createQuery(string query)
{
  OrNode *ret = new OrNode();
  int i = 0;
  int j = 0;
  i = query.find("OR");
  while(i != -1)
  {
    string str = query.substr(j, i - j);
    ret->addAnd(findAnd(str));
    i += 3;
    j = i;
    i = query.find("OR", j);
  }
  string end = query.substr(j);
  ret->addAnd(findAnd(end));
  return ret;
}

AndNode* findAnd(string query)
{
  AndNode *ret = new AndNode();
  int i = 0;
  int j = 0;
  i = query.find("AND");

  while(i != -1)
  {
    string str = query.substr(j, i - j);
    ret->addClause(findClause(str));
    i += 4;
    j = i;
    i = query.find("AND", j);
  }
  string end = query.substr(j);
  ret->addClause(findClause(end));
  return ret;
}

FieldValue* findClause(string query)
{
  string relation;
  int i = 0;
  int j = 0;
  i = query.find(" ");
  string field = query.substr(j, i - j);
  ++i;
  j = query.find(" ", i);
  relation = query.substr(i, j - i);
  ++j;
  string value = query.substr(j);
  FieldValue *ret = new FieldValue(relation, value, field);
  return ret;
}



RealPSet* UserFunc()
{
  RealPSet* ret = UserFunction();
  return ret;
}

RealPSet* Example()
{
  RealPSet* ret = ExampleFunc();
  return ret;
}

#endif


