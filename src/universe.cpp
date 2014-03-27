/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2013, The University of Texas and Bryan Marker

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



#include "transform.h"
#include <sstream>
#include <time.h>
#include "helperNodes.h"
#include "loop.h"
#include "critSect.h"

//Print out code for all generated implementations
// This takes a while for large search spaces
#define PRINTCODE

//Save generated implementations to disk
// so they can be loaded on a separate run
#define SAVETODISK


static unsigned int CURRENTSAVEVERSION = 1;
unsigned int CurrPhase = -1;

TransMap Universe::M_trans[NUMPHASES];
TransMap Universe::M_simplifiers;
TransMap Universe::M_globSimplifiers;
TransPtrMap Universe::M_transNames;
TransNameMap Universe::M_transPtrs;
unsigned int Universe::M_transCount[NUMPHASES+2];
ConsFuncMap Universe::M_consFuncMap;

Universe::Universe() {
  m_pset = NULL;
}

void Universe::Init(PSet *seed)
{
  m_pset = seed;
  m_pset->m_isTopLevel = true;
  CurrPhase = FIRSTPHASE;
  m_pset->GlobalSimplification(M_globSimplifiers, M_simplifiers);
  m_pset->BuildSizeCache();
}

void Universe::Init(string fileName)
{
  LoadFromFile(fileName);
}

Universe::~Universe()
{
  NodeVecIter iter;
  iter = m_pset->m_inTuns.begin();
  for(; iter != m_pset->m_inTuns.end(); ++iter)
    delete *iter;
  iter = m_pset->m_outTuns.begin();
  for(; iter != m_pset->m_outTuns.end(); ++iter)
    delete *iter;
  TransPtrMapIter iter2 = M_transNames.begin();
  for(; iter2 != M_transNames.end(); ++iter2) {
    delete iter2->first;
  }
  TransMapIter iter3 = M_simplifiers.begin();
  for(; iter3 != M_simplifiers.end(); ++iter3)
    delete iter3->second;
  iter3 = M_globSimplifiers.begin();
  for(; iter3 != M_globSimplifiers.end(); ++iter3)
    delete iter3->second;
  for(int i = 0; i < NUMPHASES; ++i) {
    iter3 = M_trans[i].begin();
    for(; iter3 != M_trans[i].end(); ++iter3)
      delete iter3->second;
  }
  delete m_pset;
}

bool Universe::TakeIter(unsigned int phase)
{
  bool newOne = false;

  cout << "\tStarting iteration\n";
  cout.flush();

  m_pset->BuildSizeCache();
  
  newOne = m_pset->TakeIter(M_trans[phase], M_simplifiers);

  if (newOne && !M_globSimplifiers.empty())
    if (m_pset->GlobalSimplification(M_globSimplifiers, M_simplifiers))
      m_pset->BuildSizeCache();
  
  SanityCheck();

  cout << "\tFinishing iteration\n";
  cout.flush();

  return newOne;
}

unsigned int Universe::Expand(unsigned int numIters, unsigned int phase, CullFunction cullFunc)
{
  /*
  if (phase == ROPHASE || phase == SR1PHASE) {
    m_pset->FormSets(phase);
    cout << "Formed sets\n";
  }
  */
#if DOSOPHASE
  if (phase == SOPHASE) {
    m_pset->FormSets(phase);   
  }
#endif
  
  ClearFullyExpanded();
  SanityCheck();
  if (m_pset->GlobalSimplification(M_globSimplifiers, M_simplifiers))
    m_pset->BuildSizeCache();

  unsigned int count = 0;
  double prevAlgs = TotalCount();
  bool foundNew = true;
  CurrPhase = phase;
  while ( foundNew ) {
    time_t start, end;
    time(&start);
    foundNew = TakeIter(phase);
    unsigned int total = TotalCount();
    ++count;
    if (foundNew) {
      cout << "\tSanity checking\n";
      cout.flush();
      SanityCheck();
      cout << "\tDone sanity check\n";
      cout.flush();
    }
    time(&end);
    cout << "//Done iteration " << count << " with " 
	 << total << " algorithms";
    cout << ";   increase of " << 100.0 * (total / prevAlgs - 1 )<< "%";
    cout << ";   took " << difftime(end,start) << " seconds";
    cout << endl;
    cout.flush();
    if (!foundNew) {
      time(&start);
      Prop();
#if DOSOPHASE
      //In SOPHASE, we intentionally form sets
      // to separate pack operations so the
      // same pack buffer can be used for
      // different pieces of data, sequentially.
      //We don't want to now get rid of those sets
      if (phase < SOPHASE)
#endif
	foundNew = m_pset->MergePosses(M_simplifiers, cullFunc);
      time(&end);
      if(difftime(end,start) > 10) {
	cout << "Took " << difftime(end,start) << " seconds to merge\n";
	cout.flush();
      }
      if (foundNew) {
	total = TotalCount();
	cout << total << " algorithms now" << endl;
	cout << "\tProping and sanity checking\n";
	Prop();
	SanityCheck();
	cout << "\tDone with checking\n";
	cout.flush();
      }
    }
    prevAlgs = total;
  }
  ++CurrPhase;

  cout << "Culling\n";
  cout.flush();
  time_t start, end;
  time(&start);
  Cull();
  time(&end);
  cout << "Done culling in " << difftime(end,start) << " seconds\n";
  cout.flush();

  return count;
}

void Universe::AddToMaps(Transformation *trans) 
{
  M_transNames.insert(pair<Transformation*,string>(trans,trans->GetType()));
  if (M_transPtrs.find(trans->GetType()) != M_transPtrs.end()) {
    cout << "duplicate trans name " << trans->GetType() << endl;
    cout.flush();
    throw;
  }
  M_transPtrs.insert(pair<string,Transformation*>(trans->GetType(),trans));
}

void Universe::AddTrans(const ClassType &classType, Transformation *trans, int phase)
{
  static bool hasInited = false;
  if (!hasInited) {
    for (unsigned int i = 0; i < NUMPHASES+2; ++i)
      M_transCount[i] = 0;
    hasInited = true;
  }
  if (trans->IsSingle())
    AddToMaps(trans);
  else if (trans->IsVarRef()) {
    VarTrans *var = (VarTrans*)trans;
    if (var->IsMultiRef()) {
      MultiTrans *multi = (MultiTrans*)trans;
      TransConstVecIter iter = multi->m_trans.begin();
      for(; iter != multi->m_trans.end(); ++iter)
	AddToMaps(const_cast<Transformation*>(*iter));
    }
    else
      AddToMaps(trans);
  }
  else
    throw;
  if (phase == SIMP) {
    if (trans->IsSingle()) 
      M_transCount[NUMPHASES]++;
    else 
      M_transCount[NUMPHASES] += ((MultiTrans*)trans)->NumTransformations();
    TransMapIter mapIter = M_simplifiers.find(classType);
    if (mapIter != M_simplifiers.end()) {
      mapIter->second->push_back(trans);
    }
    else {
      TransVec *vec = new TransVec;
      vec->push_back(trans);
      M_simplifiers[classType] = vec;
    }
  }
  else if (phase == GLOBSIMP) {
    if (trans->IsSingle()) 
      M_transCount[NUMPHASES+1]++;
    else 
      M_transCount[NUMPHASES+1] += ((MultiTrans*)trans)->NumTransformations();
    TransMapIter mapIter = M_globSimplifiers.find(classType);
    if (mapIter != M_globSimplifiers.end()) {
      mapIter->second->push_back(trans);
    }
    else {
      TransVec *vec = new TransVec;
      vec->push_back(trans);
      M_globSimplifiers[classType] = vec;
    }
  }
  else if (phase < 0 || phase >= NUMPHASES) 
    throw;
  else {
    if (trans->IsSingle()) 
      M_transCount[phase]++;
    else 
      M_transCount[phase] += ((MultiTrans*)trans)->NumTransformations();
    TransMapIter mapIter = M_trans[phase].find(classType);
    if (mapIter != M_trans[phase].end()) {
      mapIter->second->push_back(trans);
    }
    else {
      TransVec *vec = new TransVec;
      vec->push_back(trans);
      M_trans[phase][classType] = vec;
    }
  }
}
void Universe::SanityCheck()
{
  m_pset->SanityCheck();
}

void Universe::Prop()
{
  m_pset->ClearBeforeProp();
  m_pset->Prop();
}

void Universe::Cull()
{
  m_pset->Cull(CurrPhase);
}

void Universe::PrintAll(int algNum)
{
  time_t start,end;
  ofstream out;
  unsigned int optGraph = 0;

  //  cout << "Printing costOutput.txt\n";
  //  out.open("costOutput.txt");
  //  Print(out, COST, 0);
  //  out.close();


#ifdef MATLAB
  out.open("totalCostOutput.m");
#else
  out.open("totalCostOutput.r");
#endif
  IndStream costOut(&out, OTHERSTREAM);
  optGraph = 0;
  time(&start);
  EvalCosts(costOut, optGraph);
  out.close();
  time(&end);
  cout << "\tCost eval took " << difftime(end,start) << " seconds\n";
  cout << "Done printing\n";
  cout.flush();
  if (optGraph > 0) {
    cout << "\tOptimal graph ( " << optGraph << " ) :\n";
#if DOELEM
    IndStream optOut(&cout,ELEMSTREAM);
#elif DOSQM || DOSM
    IndStream optOut(&cout,BLISSTREAM);
#elif DOTENSORS
    IndStream optOut(&cout,TENSORSTREAM);
#endif
    Print(optOut, optGraph);
  }
  else {
    cout << "optGraph = " << optGraph << endl;
  }

#ifdef PRINTCODE
  cout << "Printing codeOutput.txt\n";
  cout.flush();
  out.open("codeOutput.txt");
#if DOELEM
    IndStream codeOut(&out,ELEMSTREAM);
#elif DOSQM || DOSM
    IndStream codeOut(&out,BLISSTREAM);
#elif DOTENSORS
    IndStream codeOut(&out,TENSORSTREAM);
#endif
  optGraph = 0;
  time(&start);
  Print(codeOut, optGraph);
  out.close();
  time(&end);
  cout << "\tTook " << difftime(end,start) << " seconds\n";

#endif

#ifdef SAVETODISK
  if (algNum != 0) {
    cout << "Saving to disk\n";
    cout.flush();
    std::stringstream str;
    str << "GraphsOut" << algNum;
    time(&start);
    SaveToFile(str.str());
    time(&end);
    cout << "\tTook " << difftime(end,start) << " seconds\n";
  }
#endif
}

void Universe::Print(IndStream &out, unsigned int &whichGraph)
{
  unsigned int graphNum = 0;
  ++graphNum;
  PossMMapIter iter = m_pset->m_posses.begin();
  for(; iter != m_pset->m_posses.end(); ++iter) {
    Poss *poss = (*iter).second;
    poss->PrintRoot(out, graphNum, whichGraph);
  }

  *out << "numAlgs = " << TotalCount() << endl;
}

void Universe::EvalCosts(IndStream &out, unsigned int &whichGraph)
{
  unsigned int optGraph, worstGraph;
  double optCost, worstCost;

#ifdef MATLAB
  *out << "transMap = containers.Map('KeyType','uint64','ValueType','char');\n";
  TransPtrMapIter iter2 = M_transNames.begin();
  for(; iter2 != M_transNames.end(); ++iter2) {
    *out << "transMap(" << (size_t)(iter2->first) << ") = '" << iter2->second << "';\n";
  }

  *out << "cost = zeros(" << TotalCount() << ",1);\n";
  *out << "refs = cell( [" << TotalCount() << " 1]);\n";
#else
  *out << "cost = array(0,dim=c(" << TotalCount() << ",1));\n";
#endif

  Prop();
  optCost = -1;
  worstCost = -1;
  unsigned int graphNum = 0;
  ++graphNum;
  PossMMapIter iter = m_pset->m_posses.begin();
  for(; iter != m_pset->m_posses.end(); ++iter) {
    Poss *poss = (*iter).second;
    poss->EvalRoot(out, graphNum, whichGraph, optGraph, optCost, worstGraph, worstCost);
  }
    
  cout << "Opt is graph " << optGraph << endl;
  cout << "\t\t " << optCost << endl;
  cout << "\tWorst is graph " << worstGraph << endl;
  cout.flush();
    
  whichGraph = optGraph;
  
  *out << "numAlgs = " << TotalCount() << endl;
}

unsigned int Universe::TotalCount() const
{
  return m_pset->TotalCount();
}

void Universe::ClearFullyExpanded()
{
  m_pset->ClearFullyExpanded();
}

void Universe::PrintStats()
{
  cout << "The universe has\n";
  for (unsigned int i = 0; i < NUMPHASES; ++i)
    cout << "\t" << M_transCount[i] << " in phase " << i << endl;
  cout << "\t" << M_transCount[NUMPHASES] << " simplifiers\n";
  cout << "\t" << M_transCount[NUMPHASES+1] << " global simplifiers\n";
}

void Universe::RegCons(ClassType type, ConstructorFunc func)
{
  if (M_consFuncMap.find(type) != M_consFuncMap.end()) {
    cout << "dup constructor for " << type << endl;
    throw;
  }
  M_consFuncMap[type] = func;
}

Node* Universe::GetBlankClassInst(ClassType type)
{
  ConsFuncMapIter iter = M_consFuncMap.find(type);
  if (iter == M_consFuncMap.end()) {
    cout << "didn't find node type " << type << endl;
    throw;
  }
  ConstructorFunc func = iter->second;
  return func();
}

void Universe::SaveToFile(string fileName) const
{
  ofstream out;
  out.open(fileName.c_str());//, ios::binary);
  Flatten(out);
  out.close();
}

void Universe::Flatten(ofstream &out) const
{
  WRITE(CURRENTSAVEVERSION);
  unsigned int tmp = M_transNames.size();
  WRITE(tmp);
  TransPtrMapConstIter iter = M_transNames.begin();
  for(; iter != M_transNames.end(); ++iter) {
    WRITE((*iter).first);
    out << (*iter).second << endl;
    //    WRITE((*iter).second);
  }
  WRITE(CurrPhase);
  Poss::FlattenStatic(out);
  Loop::FlattenStatic(out);
  WRITE(m_pset);
  bool isLoop = m_pset->IsLoop();
  WRITE(isLoop);
  bool isCrit = m_pset->IsCritSect();
  WRITE(isCrit);
  m_pset->Flatten(out);
  WRITE(END);
}

void Universe::LoadFromFile(string fileName)
{
  ifstream in;
  in.open(fileName.c_str(), ios::binary);
  Unflatten(in);
  in.close();
}

void Universe::Unflatten(ifstream &in) 
{
  unsigned int version;
  READ(version);
  if (version != CURRENTSAVEVERSION) {
    cout << version << " vs. " 
	 << CURRENTSAVEVERSION << endl;
    throw;
  }
  unsigned int numTrans;
  READ(numTrans);
  PtrMap transMap;
  for(unsigned int i = 0; i < numTrans; ++i) {
    void *old;
    READ(old);
    string name;
    getline(in, name);
    TransNameMapIter iter = M_transPtrs.find(name);
    if (iter == M_transPtrs.end()) {
      cout << "Missing transformation "
	   << name << endl;
      throw;
    }
    transMap[old] = (*iter).second;
  }
  READ(CurrPhase);
  Poss::UnflattenStatic(in);
  Loop::UnflattenStatic(in);
  PtrMap psetMap;
  PSet *oldPset;
  READ(oldPset);
  bool isLoop;
  READ(isLoop);  
  bool isCrit;
  READ(isCrit);
  if (isLoop)
    m_pset = new Loop;

  else if (isCrit)
    throw;
#if 0
    m_pset = new CritSect;
#endif
  else
    m_pset = new PSet;
  psetMap[oldPset] = m_pset;

  PtrMap possMap;
  NodeMap nodeMap;
  SaveInfo info;
  info.transMap = &transMap;
  info.possMap = &possMap;
  info.psetMap = &psetMap;
  info.nodeMap = &nodeMap;
  
  m_pset->Unflatten(in, info);
  char tmp;
  READ(tmp);
  if (tmp != END) {
    cout << "Bad end!\n";
    throw;
  }
}
