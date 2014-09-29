/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2014, The University of Texas and Bryan Marker

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
#include "loopSupport.h"
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <climits>

#if DOLLDLA

#include "LLDLAGemm.h"
#include "DLAReg.h"
#include "madd.h"
#include "vadd.h"
#include "vmmul.h"
#include "mvmul.h"
#include "vvdot.h"
#include "driverUtils.h"
#include "debug.h"
#include "LLDLAGemmTransformations.h"
#include "smmul.h"
#include "svmul.h"
#include "runtimeEvaluation.h"
#include "loopUnrolling.h"

#define DOEMPIRICALEVAL 1
#define PRINTCOSTS 1

#define DOCOMPACTLOOPUNROLLING 0
#define DO2MUTRANSFORMATIONS 1
#define DO3MUTRANSFORMATIONS 1
#define DO16MUTRANSFORMATIONS 1
#define DOLARGEMUTRANSFORMATIONS 0

#define DOPARTIALLOOPUNROLLING 1
#define PARTIALUNROLLINGSTARTCOEF 2
#define PARTIALUNROLLINGENDCOEF 16

#if DOCOMPACTLOOPUNROLLING + DOPARTIALLOOPUNROLLING > 1
do you really want to do compact unrolling and partial unrolling?
#endif

#include <sstream>

void MuNNMuGemmResults(Type precision);
double RunExample(int algNum, RealPSet* algPSet, Type precision, string opName);
RealPSet* GemvExample(Type dataType, int m, int n);
RealPSet* MVMul2Example(Type dataType, int m, int n, int p);
RealPSet* MAdd2Example(Type dataType, int m, int n);
RealPSet* VAdd2Example(Type dataType, int m);
RealPSet* VAddExample(Type dataType, int m);
RealPSet* VMVMulExample(Type dataType, int m, int n);
RealPSet* SMMulExample(Type dataType, int m, int n);
RealPSet* VMMulExample(Type dataType, int m, int n);
RealPSet* SVMulRowExample(Type dataType, int m);
RealPSet* SVMulColExample(Type dataType, int m);
RealPSet* MVMulExample(Type dataType, int m, int n);
RealPSet* MAddExample(Type dataType, int m, int n);
RealPSet* DotExample(Type dataType, int m);
RealPSet* GemmExample(Type dataType, int m, int n, int p);
RealPSet* DoubleGemmExample(Type dataType, int m, int n, int p, int k);

Trans transA, transB;

Architecture* arch;

ImplementationMap* ImpStrMap(Universe *uni)
{
  ImplementationMap* impMap = new ImplementationMap();
  GraphNum i;
  for (i = 1; i <= uni->TotalCount(); i++) {
    std::stringbuf sbuf;
    std::ostream out(&sbuf);
    IndStream istream = IndStream(&out, LLDLASTREAM);
    uni->Print(istream, i);
    impMap->insert(NumImplementationPair(i, sbuf.str()));
  }
  return impMap;
}

void PrintImpMap(ImplementationRuntimeMap &impTimes)
{
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    cout << "IMPLEMENTATION # " << std::to_string((long long int) mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      cout << std::to_string((long double) *vit) << endl;
    }
    cout << endl;
  }
}

double BestFlopsPerCycle(Type type, ImplementationRuntimeMap &impTimes, double flopCost) {
  double peakFlopsPerCycle = arch->FlopsPerCycle(type);
  double bestFlopsPerCycle = 0;
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      double totalTimeInCycles = *vit;
      double actualFlopsPerCycle = flopCost / totalTimeInCycles;
      double pctPeak = (actualFlopsPerCycle / peakFlopsPerCycle) * 100;
      if (actualFlopsPerCycle > bestFlopsPerCycle) {
	bestFlopsPerCycle = actualFlopsPerCycle;
      }
      if (pctPeak > 100) {
	cout << "pctPeak > 100\n";
	throw;
      }
    }
  }
  return bestFlopsPerCycle;
}

GraphNum PrintImpMapInFlops(Type type, ImplementationRuntimeMap &impTimes, double flopCost) {
  double peakFlopsPerCycle = arch->FlopsPerCycle(type);
  GraphNum bestImpNum = 0;
  double bestFlopsPerCycle = 0;
  ImplementationRuntimeMapIter mit;
  for (mit = impTimes.begin(); mit != impTimes.end(); ++mit) {
    TimeVecIter vit;
    cout << "IMPLEMENTATION # " << std::to_string((long long int) mit->first) << endl;
    for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
      double totalTimeInCycles = *vit;
      double actualFlopsPerCycle = flopCost / totalTimeInCycles;
      double pctPeak = (actualFlopsPerCycle / peakFlopsPerCycle) * 100;
      if (actualFlopsPerCycle > bestFlopsPerCycle) {
	bestFlopsPerCycle = actualFlopsPerCycle;
	bestImpNum = mit->first;
      }
      cout << "Flops per cycle = " << std::to_string((long double) actualFlopsPerCycle);
      cout << "\t%Peak = " << std::to_string((long double) pctPeak) << endl;
      if (pctPeak > 100) {
	cout << "pctPeak > 100\n";
	throw;
      }
    }
    cout << endl;
  }
  cout << "Best flops/cycle achieved: " << std::to_string((long double) bestFlopsPerCycle / 1.0e9) << endl;
  cout << "Best percent of peak: " << std::to_string((long double) (bestFlopsPerCycle / peakFlopsPerCycle) * 100) << endl;
  return bestImpNum;      
}

void AddGemmTrans()
{
    // Convert gemm into loop over mvmul
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToMVMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  // Transform gemm into loop over vmmuls
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmToVMMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  //Introduces loops in the m, n, and k dimensions, respectively
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

#if DOCOMPACTLOOPUNROLLING
#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
#endif

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
#endif
#endif // DOCOMPACTLOOPUNROLLING


#if DOCOMPACTLOOPUNROLLING
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);

#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
#endif

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new LLDLAGemmLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
#endif

#endif // DOCOMPACTLOOPUNROLLING

  return;
}

void AddVVDotTrans()
{
  Universe::AddTrans(VVDot::GetClass(), new VVDotToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddMAddTrans()
{
  Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuDouble), LLDLALOOPPHASE);

  Universe::AddTrans(MAdd::GetClass(), new MAddLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble), LLDLALOOPPHASE);

  Universe::AddTrans(MAdd::GetClass(), new MAddToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddMVMulTrans()
{
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuDouble), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddSMMulTrans()
{
  //Introduces loops in the m and n dimension for SMMul
  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(SMMul::GetClass(), new SMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble), LLDLALOOPPHASE);

  return;
}

void AddUnrollingTrans()
{

#if DOCOMPACTLOOPUNROLLING
#if DO2MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(),
		     new CompactlyUnrollLoop(2), LLDLALOOPUNROLLPHASE);
#endif // DO2MUTRANSFORMATIONS

#if DO3MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
		     new CompactlyUnrollLoop(3), LLDLALOOPUNROLLPHASE);
#endif // DO3MUTRANSFORMATIONS

#if DO16MUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
		     new CompactlyUnrollLoop(16), LLDLALOOPUNROLLPHASE);
#endif // DO3MUTRANSFORMATIONS

#if DOLARGEMUTRANSFORMATIONS
  Universe::AddTrans(SplitSingleIter::GetClass(), 
		     new CompactlyUnrollLoop(bigSize / LLDLA_MU), LLDLALOOPUNROLLPHASE);
#endif // DOLARGEMUTRANSFORMATIONS

#endif // DOCOMPACTLOOPUNROLLING

#if DOPARTIALLOOPUNROLLING
  for (unsigned int mult = PARTIALUNROLLINGSTARTCOEF; mult <= PARTIALUNROLLINGENDCOEF; mult += 2) {
    Universe::AddTrans(SplitSingleIter::GetClass(), new PartiallyUnrollLoop(mult), LLDLALOOPUNROLLPHASE);
  }
#endif // DOPARTIALLOOPUNROLLING

}

void AddSVMulTrans()
{
  Universe::AddTrans(SVMul::GetClass(), new SVMulToRegArith(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);

  Universe::AddTrans(SVMul::GetClass(), new SVMulToRegArith(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);

  return;
}

void AddVMMulTrans()
{
  // Transformers for vector matrix multiply
  Universe::AddTrans(VMMul::GetClass(), new VMMulToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMK, LLDLAMuDouble), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble), LLDLALOOPPHASE);

  return;
}

void AddVAddTrans()
{
  Universe::AddTrans(VAdd::GetClass(), new VAddToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddTrans()
{
  AddGemmTrans();
  AddVVDotTrans();
  AddMAddTrans();
  AddMVMulTrans();
  AddSMMulTrans();
  AddSVMulTrans();
  AddVMMulTrans();
  AddVAddTrans();

  AddUnrollingTrans();
  
}

void Usage()
{
  cout << "./driver arg1 arg2 ...\n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Gemm  N/T N/T F/D M N P\n";
  cout <<"         2  -> Double Gemm  N/T N/T F/D M N P K\n";
  cout <<"         3  -> Dot prod F/D M\n";
  cout <<"         4  -> Matrix add F/D M N\n";
  cout <<"         5  -> Matrix vector multiply F/D M N\n";
  cout <<"         6  -> Scalar column vector multiply F/D M\n";
  cout <<"         7  -> Scalar row vector multiply F/D M\n";
  cout <<"         8  -> Vector matrix multiply F/D M N\n";
  cout <<"         9  -> Scalar matrix multiply F/D M N\n";
  cout <<"        10  -> Vector add F/D M\n";
  cout <<"        11  -> Vector add twice F/D M\n";
  cout <<"        12  -> Vector matrix vector multiply F/D M N\n";
  cout <<"        13  -> Matrix add twice F/D M N\n";
  cout <<"        14  -> Matrix vector multiply twice F/D M N P\n";
  cout <<"        15  -> Gemv F/D M N\n";
  cout <<"        16  -> Run (mu x n) (n x mu) tests with precision F/D\n";
}

int main(int argc, const char* argv[])
{
#ifdef _OPENMP
  omp_set_num_threads(1);
  omp_set_nested(true);
#endif

  int m, n, p, k;
  Type precision;

  arch = new HaswellMacbook();


  //  PrintType printType = CODE;

  RealPSet* algPSet;
  //  GraphNum whichGraph = 0;
  int algNum;
  string opName;


  if(argc < 2) {
    Usage();
    return 0;
  }
  else {
    algNum = atoi(argv[1]);
    switch(algNum) {
    case(1):
      if (argc != 8) {
	Usage();
	return 0;
      }
      opName = "dxt_gemm";
      transA = CharToTrans(*argv[2]);
      transB = CharToTrans(*argv[3]);
      precision = CharToType(*argv[4]);
      m = atoi(argv[5]);
      n = atoi(argv[6]);
      p = atoi(argv[7]);
      algPSet = GemmExample(precision, m, n, p);
      break;
    case(2):
      if (argc != 9) {
	Usage();
	return 0;
      }
      opName = "dxt_double_gemm";
      precision = CharToType(*argv[4]);
      m = atoi(argv[5]);
      n = atoi(argv[6]);
      p = atoi(argv[7]);
      k = atoi(argv[8]);
      algPSet = DoubleGemmExample(precision, m, n, p, k);
      transA = CharToTrans(*argv[2]);
      transB = CharToTrans(*argv[3]);
      break;
    case(3):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_dot";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      algPSet = DotExample(precision, m);
      break;
    case(4):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_madd";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = MAddExample(precision, m, n);
      break;
    case(5):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_mvmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = MVMulExample(precision, m, n);
      break;
    case(6):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_sv_col_mul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      algPSet = SVMulColExample(precision, m);
      break;
    case(7):
      if (argc != 4) {
	Usage();
	return 0;
      }
      precision = CharToType(*argv[2]);
      m = atoi(argv[4]);
      opName = "dxt_sv_row_mul";
      algPSet = SVMulRowExample(precision, m);
      break;
    case(8):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_vmmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = VMMulExample(precision, m, n);
      break;
    case(9):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_smmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = SMMulExample(precision, m, n);
      break;
    case(10):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_vadd";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      algPSet = VAddExample(precision, m);
      break;
    case(11):
      if (argc != 4) {
	Usage();
	return 0;
      }
      opName = "dxt_vadd2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      algPSet = VAdd2Example(precision, m);
      break;
    case(12):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_vmvmul";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = VMVMulExample(precision, m, n);
      break;
    case(13):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_madd2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = MAdd2Example(precision, m, n);
      break;
    case(14):
      if (argc != 6) {
	Usage();
	return 0;
      }
      opName = "dxt_mvmul2";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      p = atoi(argv[5]);
      algPSet = MVMul2Example(precision, m, n, p);
      break;
    case(15):
      if (argc != 5) {
	Usage();
	return 0;
      }
      opName = "dxt_gemv";
      precision = CharToType(*argv[2]);
      m = atoi(argv[3]);
      n = atoi(argv[4]);
      algPSet = GemvExample(precision, m, n);
      break;
    case(16):
      if (argc != 3) {
	Usage();
	return 0;
      }
      precision = CharToType(*argv[2]);
      MuNNMuGemmResults(precision);
      return 0;
    default:
      Usage();
      return 0;
    }
  }

  RunExample(algNum, algPSet, precision, opName);
  return 0;
}

void MuNNMuGemmResults(Type precision) {
  int m = arch->VecRegWidth(precision);
  int n = m;
  int p = 16;
  int p_inc = 16;
  int num_trials = 15;
  int algNum = 1;
  string opName = "dxt_gemm";
  string resultFileName = "dxt_gemm_mu_by_n_times_n_by_mu_trials.csv";
  std::ofstream outFile(resultFileName);

  outFile << "m, n, p, bestFlopsPerCylcle\n";

  for (int i = 0; i < num_trials; i++) {
    RealPSet* algPSet = GemmExample(precision, m, n, p);
    double bestFlops = RunExample(algNum, algPSet, precision, opName);
    outFile << std::to_string((long long int) m) << ", ";
    outFile << std::to_string((long long int) n) << ", ";
    outFile << std::to_string((long long int) p) << ", ";
    outFile << std::to_string((double) bestFlops) << "\n";
    p += p_inc;
  }
  outFile.close();
  return;
}

double RunExample(int algNum, RealPSet* algPSet, Type precision, string opName)
{
  RegAllLLDLANodes();
  AddTrans();

  int numIters = -1;
  Cost flopCost = 0;
  Universe uni;
  time_t start, start2, end;
  string absImpStr;

  uni.PrintStats();

  cout << "Creating startSet\n";

  RealPSet *startSet = algPSet;
  
  cout << "Created startSet\n";

  uni.Init(startSet);
  
  cout << "Initialized universe\n";
  
  uni.Prop();
  GraphIter* graphIter = new GraphIter(startSet->m_posses.begin()->second);
  cout << "Printing evaluation code\n";
  flopCost = graphIter->EvalAndSetBest();
  // Print abstract implementation to string for use in testing
  // EXTREMELY HACKY, I could not figure out how to redirect an
  // ostream to a string
  std::stringstream ss;
  IndStream optOut(&ss, LLDLASTREAM);
  graphIter->PrintRoot(optOut, 0, true, startSet);
  absImpStr = ss.str();

  cout << "IMPLEMENTATION FOR CORRECTNESS CHECK:\n" << absImpStr;
  cout << "Flops for operation = " << std::to_string((long double) flopCost) << endl;

  time(&start);

#if DOLLDLALOOPPHASE
  if (CurrPhase == LLDLALOOPPHASE) {

    cout << "Expanding LLDLA loop phase\n";

    uni.Expand(-1, LLDLALOOPPHASE, LLDLACull);
    time(&end);

    cout << "LLDLALOOP phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";

    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);

    cout << "Propagation took " << difftime(end,start2) << " seconds\n";

  }
#endif


#if DOLLDLALOOPUNROLLPHASE
  if (CurrPhase == LLDLALOOPUNROLLPHASE) {
    cout << "LLDLALOOPUNROLL phase\n";
    uni.Expand(-1, LLDLALOOPUNROLLPHASE, LLDLACull);
    time(&end);
    cout << "LLDLALOOPUNROLL phase took " << difftime(end,start) << " seconds\n";
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOLLDLAPRIMPHASE
  if (CurrPhase == LLDLAPRIMPHASE) {
    cout << "Expanding LL DLA prim phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, LLDLAPRIMPHASE, LLDLACull);
    time(&end);
    cout << "LLDLAPRIM phase took " << difftime(end,start2) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

  cout << "Full expansion took " << difftime(end,start) << " seconds\n";
  cout.flush();

#if DOEMPIRICALEVAL  
  cout << "Writing all implementations to runtime eval files\n";

  int numIterations = 10;
  RuntimeTest rtest(precision, opName, uni.m_argNames, uni.m_declarationVectors, uni.m_constantDefines, numIterations);
  string evalDirName = "runtimeEvaluation";
  RuntimeEvaluator evaler = RuntimeEvaluator(evalDirName);
  cout << "About to evaluate\n";
  ImplementationRuntimeMap impMap = evaler.EvaluateImplementationsWithCorrectnessCheck(rtest, ImpStrMap(&uni), absImpStr);


  cout << "Done evaluating\n";
  GraphNum best = PrintImpMapInFlops(precision, impMap, flopCost);
  cout << "All implementations printed\n";
  cout << "Best times";
  
#endif //DOEMPIRICALEVAL

#if 1
  uni.PrintAll(algNum, best);
#else
  uni.PrintBest();
#endif

#if PRINTCOSTS
  uni.PrintCosts(impMap);
#endif

  double bestFPS = BestFlopsPerCycle(precision, impMap, flopCost);
  return bestFPS;
}

RealPSet* GemvExample(Type dataType, int m, int n)
{
  InputNode* xIn = new InputNode("x input", n, 1, "X",
				 1, n,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);

  InputNode* AIn = new InputNode("a input", m, n, "A",
				 1, m,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* alphaIn = new InputNode("alpha input", 1, 1, "Alpha",
				     1, m,
				     "AlphaNumRows", "AlphaNumCols",
				     "AlphaRowStride", "AlphaColStride", dataType);

  InputNode* betaIn = new InputNode("beta input", 1, 1, "Beta",
				    1, m,
				    "BetaNumRows", "BetaNumCols",
				    "BetaRowStride", "BetaColStride", dataType);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  Tunnel* tunAlpha = new Tunnel(POSSTUNIN);
  tunAlpha->AddInput(alphaIn, 0);

  Tunnel* tunBeta = new Tunnel(POSSTUNIN);
  tunBeta->AddInput(betaIn, 0);

  SVMul* by = new SVMul(COLVECTOR, ABSLAYER);
  by->AddInputs(4,
		tunBeta, 0,
		tunY, 0);

  MVMul* axMul = new MVMul(ABSLAYER);
  axMul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunZ, 0);

  SVMul* alphaAXMul = new SVMul(COLVECTOR, ABSLAYER);
  alphaAXMul->AddInputs(4,
			tunAlpha, 0,
			axMul, 0);

  VAdd* sumVecs = new VAdd(COLVECTOR, ABSLAYER);
  sumVecs->AddInputs(4,
		     alphaAXMul, 0,
		     by, 0);

  Poss* innerPoss = new Poss(sumVecs, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MVMul2Example(Type dataType, int m, int n, int p)
{
  InputNode* xIn = new InputNode("x input", p, 1, "X",
				 1, p,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", n, 1, "Y",
				 1, n,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);

  InputNode* AIn = new InputNode("a input", m, n, "A",
				 1, m,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* BIn = new InputNode("b input", n, p, "B",
				 1, n,
				 "BNumRows", "BNumCols",
				 "BRowStride", "BColStride", dataType);
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(AIn, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(BIn, 0);

  MVMul* mvmul1 = new MVMul(ABSLAYER);
  mvmul1->AddInputs(6,
		    tunB, 0,
		    tunX, 0,
		    tunY, 0);

  MVMul* mvmul2 = new MVMul(ABSLAYER);
  mvmul2->AddInputs(6,
		    tunA, 0,
		    mvmul1, 0,
		    tunZ, 0);

  Poss* innerPoss = new Poss(mvmul2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MAdd2Example(Type dataType, int m, int n)
{
  InputNode* xIn = new InputNode("x input", m, n, "X",
				 1, m,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", m, n, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* zIn = new InputNode("z input", m, n, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  MAdd* madd1 = new MAdd(ABSLAYER);
  madd1->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  MAdd* madd2 = new MAdd(ABSLAYER);
  madd2->AddInputs(4,
		   tunZ, 0,
		   madd1, 0);

  Poss* innerPoss = new Poss(madd2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VAdd2Example(Type dataType, int m)
{
  InputNode* xIn = new InputNode("x input", m, 1, "X",
				 1, m,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  VAdd* vadd1 = new VAdd(COLVECTOR, ABSLAYER);
  vadd1->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  VAdd* vadd2 = new VAdd(COLVECTOR, ABSLAYER);
  vadd2->AddInputs(4,
		   tunZ, 0,
		   vadd1, 0);

  Poss* innerPoss = new Poss(vadd2, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

}

RealPSet* VAddExample(Type dataType, int m)
{
  InputNode* xIn = new InputNode("x input", m, 1, "X",
				 1, m,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);
  
  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  VAdd* vadd = new VAdd(COLVECTOR, ABSLAYER);
  vadd->AddInputs(4,
		  tunX, 0,
		  tunY, 0);

  Poss* innerPoss = new Poss(vadd, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;

}

RealPSet* VMVMulExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A",
				 1, m,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* xIn = new InputNode("x input", n, 1, "X",
				 1, n,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* zIn = new InputNode("z input", m, 1, "Z",
				 1, m,
				 "ZNumRows", "ZNumCols",
				 "ZRowStride", "ZColStride", dataType);

  InputNode* yIn = new InputNode("y input", 1, m, "Y",
				 1, 1,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  InputNode* wIn = new InputNode("w input", 1, 1, "W",
				 1, 1,
				 "WNumRows", "WNumCols",
				 "WRowStride", "WColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunZ = new Tunnel(POSSTUNIN);
  tunZ->AddInput(zIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  Tunnel* tunW = new Tunnel(POSSTUNIN);
  tunW->AddInput(wIn, 0);

  MVMul* mvmul = new MVMul(ABSLAYER);
  mvmul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunZ, 0);

  VVDot* vvdot = new VVDot(ABSLAYER);
  vvdot->AddInputs(6,
		   tunY, 0,
		   mvmul, 0,
		   tunW, 0);

  Poss* innerPoss = new Poss(vvdot, true);
  RealPSet* innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SMMulExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A",
				 n, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);
  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 1, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SMMul* smmul = new SMMul(ABSLAYER);
  smmul->AddInputs(4,
		   tunX, 0,
		   tunA, 0);

  Poss *innerPoss = new Poss(smmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* VMMulExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A",
				 n, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);
  InputNode* xIn = new InputNode("x input", 1, m, "X",
				 m, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  InputNode* yIn = new InputNode("y input", 1, n, "Y",
				 n, 1,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  VMMul* vmmul = new VMMul(ABSLAYER);
  vmmul->AddInputs(6,
		   tunX, 0,
		   tunA, 0,
		   tunY, 0);

  Poss *innerPoss = new Poss(vmmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SVMulRowExample(Type dataType, int m)
{
  InputNode* Ain = new InputNode("A input", 1, m, "A",
				 m, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 m, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SVMul* svmul = new SVMul(ROWVECTOR, ABSLAYER);
  svmul->AddInputs(4,
		   tunX, 0,
		   tunA, 0);

  Poss *innerPoss = new Poss(svmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* SVMulColExample(Type dataType, int m)
{
  InputNode* Ain = new InputNode("A input", m, 1, "A",
				 m, 1,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* xIn = new InputNode("x input", 1, 1, "X",
				 m, 1,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  SVMul* svmul = new SVMul(COLVECTOR, ABSLAYER);
  svmul->AddInputs(4,
		   tunX, 0,
		   tunA, 0);

  Poss *innerPoss = new Poss(svmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MVMulExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A",
				 1, m,
				 "ANumRows", "ANumCols",
				 "ARowStride", "AColStride", dataType);

  InputNode* xIn = new InputNode("x input", n, 1, "X",
				 1, n,
				 "XNumRows", "XNumCols",
				 "XRowStride", "XColStride", dataType);
  InputNode* yIn = new InputNode("y input", m, 1, "Y",
				 1, m,
				 "YNumRows", "YNumCols",
				 "YRowStride", "YColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunX = new Tunnel(POSSTUNIN);
  tunX->AddInput(xIn, 0);

  Tunnel* tunY = new Tunnel(POSSTUNIN);
  tunY->AddInput(yIn, 0);

  MVMul* mvmul = new MVMul(ABSLAYER);
  mvmul->AddInputs(6,
		   tunA, 0,
		   tunX, 0,
		   tunY, 0);

  Poss *innerPoss = new Poss(mvmul, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* MAddExample(Type dataType, int m, int n)
{
  InputNode* Ain = new InputNode("A input", m, n, "A", 
				 1, m,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride", dataType);

  InputNode* Bin = new InputNode("B input", m, n, "B", 
				 1, m,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride", dataType);

  Tunnel* tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel* tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin, 0);

  MAdd* madd = new MAdd(ABSLAYER);
  madd->AddInputs(4,
		 tunA, 0,
		 tunB, 0);

  Poss *innerPoss = new Poss(madd, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* DotExample(Type dataType, int m)
{
  InputNode* Ain = new InputNode("A input", 1, m, "A", 
				 m, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride", dataType);

  InputNode* Bin = new InputNode("B input", m, 1, "B", 
				 m, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride", dataType);

  InputNode* Cin = new InputNode("C input", 1, 1, "C", 
				 m, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride", dataType);

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  VVDot* dot = new VVDot(ABSLAYER);
  dot->AddInputs(6,
		 tunA, 0,
		 tunB, 0,
		 tunC, 0);

  Poss *innerPoss = new Poss(dot, true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0), 0);

  Poss *outerPoss = new Poss(Cout, true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* GemmExample(Type dataType, int m, int n, int p)
{
  InputNode *Ain= new InputNode("A input", m, p, "A",
				 p, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride", dataType);

  InputNode *Bin = new InputNode("B input", p, n, "B",
				 n, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride", dataType);

  InputNode *Cin = new InputNode("C input", m, n, "C",
				 n, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride", dataType);

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain, 0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin, 0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin, 0);

  Gemm *gemm = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, dataType);
  gemm->AddInputs(6,
		  tunA, 0,
		  tunB, 0,
		  tunC, 0);

  Poss *innerPoss = new Poss(gemm,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}

RealPSet* DoubleGemmExample(Type dataType, int m, int n, int p, int k)
{
  InputNode *Ain = new InputNode("A input",  m, n, "A",
				 n, 1,
				 "ANumRows","ANumCols",
				 "ARowStride","AColStride", dataType);

  InputNode *Bin = new InputNode("B input", n, p, "B",
				 p, 1,
				 "BNumRows","BNumCols",
				 "BRowStride","BColStride", dataType);

  InputNode *Cin = new InputNode("C input",  p, k, "C",
				 k, 1,
				 "CNumRows","CNumCols",
				 "CRowStride","CColStride", dataType);

  Tunnel *tunA = new Tunnel(POSSTUNIN);
  tunA->AddInput(Ain,0);

  Tunnel *tunB = new Tunnel(POSSTUNIN);
  tunB->AddInput(Bin,0);

  Tunnel *tunC = new Tunnel(POSSTUNIN);
  tunC->AddInput(Cin,0);

  Gemm *gemm1 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, dataType);
  gemm1->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  tunC,0);

  Gemm *gemm2 = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, dataType);
  gemm2->AddInputs(6,
		  tunA,0,
		  tunB,0,
		  gemm1,0);

  Poss *innerPoss = new Poss(gemm2,true);
  RealPSet *innerSet = new RealPSet(innerPoss);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Cout,true);
  RealPSet *outerSet = new RealPSet(outerPoss);
  
  return outerSet;
}


#endif //DOLLDLA
