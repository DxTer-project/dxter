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


#if DORQO

RealPSet* Example1();
RealPSet* Example2();
RealPSet* Example3();
RealPSet* Example4();
RealPSet* Example5();
void Example1Run();

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
  
#if 0
  uni.PrintAll(algNum);
#else
  uni.PrintBest();
#endif


  LOG_END();

  Example1Run();
  return 0;
}

RealPSet* Example1()
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

  HJoin *join = new HJoin("u", joinFields0, joinFields1);

  join->AddInput(inA, 0);
  join->AddInput(inB, 0);


  vector<string> joinFields2;
  joinFields2.push_back("y");

  vector<string> joinFields3;
  joinFields3.push_back("c");

  HJoin *join2 = new HJoin("b", joinFields2, joinFields3);

  join2->AddInput(join, 0);
  join2->AddInput(inC, 0);

  set<string> projFields;
  projFields.insert("x");
  projFields.insert("u");
  projFields.insert("b");

  Projection *proj = new Projection("x", projFields);

  proj->AddInput(join2, 0);

  Poss *poss = new Poss(1, proj);
  RealPSet *pset = new RealPSet(poss);

  Example1Run();
  return pset;
}

void Example1Run()
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
}

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


#endif


