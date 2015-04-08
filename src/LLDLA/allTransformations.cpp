#include "allTransformations.h"

#if DOLLDLA

#include "copy.h"
#include "copyColLoopRef.h"
#include "copyRowLoopRef.h"
#include "copyToContigCopy.h"
#include "eliminateRecombine.h"
#include "eliminateRecombinePartition.h"
#include "eliminateStoreLoad.h"
#include "LLDLATranspose.h"
#include "loadToContigLoad.h"
#include "mmul.h"
#include "mmulTransformations.h"
#include "madd.h"
#include "mvmul.h"
#include "mvmulToSVMul.h"
#include "mvmulToVVDot.h"
#include "mvmulPack.h"
#include "pack.h"
#include "packToCopyAndZero.h"
#include "partition.h"
#include "recombine.h"
#include "regLoadStore.h"
#include "residualSVMulToRegArith.h"
#include "residualVAddToRegArith.h"
#include "residualVVDotToRegArith.h"
#include "setToZero.h"
#include "setToZeroColLoopRef.h"
#include "setToZeroRowLoopRef.h"
#include "smmul.h"
#include "svmul.h"
#include "svmulPackResidualToVRW.h"
#include "svmulSplitToMainAndResidual.h"
#include "unpack.h"
#include "unpackToPartAndCopy.h"
#include "vmmul.h"
#include "vadd.h"
#include "vaddPackResidualToVRW.h"
#include "vaddSplitToMainAndResidual.h"
#include "vrwSVMulToRegArith.h"
#include "vrwVAddToRegArith.h"
#include "vrwVVDotToRegArith.h"
#include "vvdot.h"
#include "vvdotPackToMultipleOfMu.h"
#include "vvdotSplitToMainAndResidual.h"

#define DOCOMPACTLOOPUNROLLING 0
#define DO2MUTRANSFORMATIONS 0
#define DO3MUTRANSFORMATIONS 0
#define DO16MUTRANSFORMATIONS 0
#define DOLARGEMUTRANSFORMATIONS 0

#define DOPARTIALLOOPUNROLLING 0
#define PARTIALUNROLLINGSTARTCOEF 0
#define PARTIALUNROLLINGENDCOEF 0

#if DOCOMPACTLOOPUNROLLING + DOPARTIALLOOPUNROLLING > 1
do you really want to do compact unrolling and partial unrolling?
#endif


void AddGemmTrans() {
  // Convert gemm into loop over mvmul
  Universe::AddTrans(Gemm::GetClass(), new MMulToMVMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  // Transform gemm into loop over vmmuls
  Universe::AddTrans(Gemm::GetClass(), new MMulToVMMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  //Introduces loops in the m, n, and k dimensions, respectively
  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

  if (DOCOMPACTLOOPUNROLLING) {
    if (DO2MUTRANSFORMATIONS) {
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA2MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
    }

    if (DO3MUTRANSFORMATIONS) {
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA3MuSingle, REAL_SINGLE), LLDLALOOPPHASE);
    }
  }


    if (DOCOMPACTLOOPUNROLLING) {
    Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
    Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
    Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle, REAL_DOUBLE), LLDLALOOPPHASE);

    if (DO2MUTRANSFORMATIONS) {
    Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
    Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
    Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA2MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
  }

    if(DO3MUTRANSFORMATIONS) {
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
      Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLA3MuSingle, REAL_DOUBLE), LLDLALOOPPHASE);
    }

  }

    return;
  }

void AddVVDotTrans() {
  Universe::AddTrans(VVDot::GetClass(), new VVDotSplitToMainAndResidual(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(VVDot::GetClass(), new VVDotToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(VVDot::GetClass(), new ResidualVVDotToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(VVDot::GetClass(), new VRWVVDotToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);


  return;
}

void AddMAddTrans() {
    Universe::AddTrans(MAdd::GetClass(), new MAddToVAddLoopRef(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);
    Universe::AddTrans(MAdd::GetClass(), new MAddToVAddLoopRef(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);

    return;
}

void AddMVMulTrans() {
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuDouble), LLDLALOOPPHASE);
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble), LLDLALOOPPHASE);
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle), LLDLALOOPPHASE);
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulToSVMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulToVVDot(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  //Universe::AddTrans(MVMul::GetClass(), new MVMulToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  //  Universe::AddTrans(MVMul::GetClass(), new MVMulPackOutput(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddSMMulTrans() {
  Universe::AddTrans(SMMul::GetClass(), new SMulToSVMul(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);
  Universe::AddTrans(SMMul::GetClass(), new SMulToSVMul(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);

  return;
}

void AddUnrollingTrans() {

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

void AddSVMulTrans() {
    Universe::AddTrans(SVMul::GetClass(), new SVMulLoopRef(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);
    Universe::AddTrans(SVMul::GetClass(), new SVMulLoopRef(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);

    Universe::AddTrans(SVMul::GetClass(), new VRWSVMulToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

    Universe::AddTrans(SVMul::GetClass(), new ResidualSVMulToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

    Universe::AddTrans(SVMul::GetClass(), new SVMulSplitToMainAndResidual(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);
    Universe::AddTrans(SVMul::GetClass(), new SVMulSplitToMainAndResidual(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);

    return;
}

void AddVMMulTrans() {
  Universe::AddTrans(VMMul::GetClass(), new VMMulToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMK, LLDLAMuDouble), LLDLALOOPPHASE);

  Universe::AddTrans(VMMul::GetClass(), new VMMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble), LLDLALOOPPHASE);

  return;
}

void AddVAddTrans() {
  Universe::AddTrans(VAdd::GetClass(), new VAddLoopRef(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);
  Universe::AddTrans(VAdd::GetClass(), new VAddLoopRef(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);

  Universe::AddTrans(VAdd::GetClass(), new VAddSplitToMainAndResidual(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);
Universe::AddTrans(VAdd::GetClass(), new VAddSplitToMainAndResidual(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);

Universe::AddTrans(VAdd::GetClass(), new ResidualVAddToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(VAdd::GetClass(), new VRWVAddToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddTransposeTrans() {
  Universe::AddTrans(LLDLATranspose::GetClass(), new LLDLATransposeLowerLayer(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);
  Universe::AddTrans(LLDLATranspose::GetClass(), new LLDLATransposeLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER), LLDLALOOPPHASE);

  return;
}

void AddPartitionRecombineTrans() {
  Universe::AddTrans(Partition::GetClass(), new PartitionLowerLayer(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);
  Universe::AddTrans(Partition::GetClass(), new PartitionLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER), LLDLALOOPPHASE);

  Universe::AddTrans(Recombine::GetClass(), new RecombineLowerLayer(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);
  Universe::AddTrans(Recombine::GetClass(), new RecombineLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER), LLDLALOOPPHASE);
}

void AddSetToZeroTrans() {
  Universe::AddTrans(SetToZero::GetClass(), new SetToZeroLowerLayer(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(SetToZero::GetClass(), new SetToZeroLowerLayer(LLDLAMIDLAYER, LLDLAPRIMITIVELAYER), LLDLALOOPPHASE);

  Universe::AddTrans(SetToZero::GetClass(), new SetToZeroRowLoopRef(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(SetToZero::GetClass(), new SetToZeroColLoopRef(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);

}

void AddPackTrans() {
  Universe::AddTrans(Pack::GetClass(), new PackToCopyAndZero(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);
}

void AddUnpackTrans() {
  Universe::AddTrans(Unpack::GetClass(), new UnpackToPartAndCopy(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);
}

void AddCopyTrans() {
  Universe::AddTrans(Copy::GetClass(), new CopyToContigCopy(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(Copy::GetClass(), new CopyColLoopRef(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(Copy::GetClass(), new CopyRowLoopRef(ABSLAYER, LLDLAMIDLAYER), LLDLALOOPPHASE);
}

void AddRTLOptimizations() {
  Universe::AddTrans(StoreFromRegs::GetClass(), new EliminateStoreLoad(ABSLAYER, ABSLAYER), SIMP);
  Universe::AddTrans(Recombine::GetClass(), new EliminateRecombinePartition(ABSLAYER, ABSLAYER), SIMP);
  Universe::AddTrans(Recombine::GetClass(), new EliminateRecombine(ABSLAYER, ABSLAYER), SIMP);
  Universe::AddTrans(LoopTunnel::GetClass(), new HoistLoad, SIMP);
}

void AddPrimPhaseConversions() {
  Universe::AddTrans(LoadToRegs::GetClass(), new LoadToContigLoad(ABSLAYER, ABSLAYER), LLDLAPRIMPHASE);
}

void AddTransformations() {
  AddGemmTrans();
  AddVVDotTrans();
  AddMAddTrans();
  AddMVMulTrans();
  AddSMMulTrans();
  AddSVMulTrans();
  AddVMMulTrans();
  AddVAddTrans();

  AddRTLOptimizations();
  AddPrimPhaseConversions();

  AddTransposeTrans();
  AddUnrollingTrans();
}

#endif // DOLLDLA
