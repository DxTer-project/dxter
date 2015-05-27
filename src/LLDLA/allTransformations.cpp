#include "allTransformations.h"

#if DOLLDLA

#include "addMulToFMA.h"
#include "copy.h"
#include "copyColLoopRef.h"
#include "copyRowLoopRef.h"
#include "copyToContigCopy.h"
#include "distributeMVMulOverMAdd.h"
#include "eliminatePackedStore.h"
#include "eliminatePackedStoreLoad.h"
#include "eliminateRecombine.h"
#include "eliminateRecombinePartition.h"
#include "eliminateStore.h"
#include "eliminateStoreLoad.h"
#include "eliminateStoreLoadOut.h"
#include "hoistDuplicateLoad.h"
#include "hoistLoadToRegs.h"
#include "LLDLATranspose.h"
#include "loadToContigLoad.h"
#include "mmul.h"
#include "mmulSplitAlongP.h"
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
#include "svmulAdd.h"
#include "svmulAddLoopRef.h"
#include "svmulAddSplitToMainAndResidual.h"
#include "residualSVMulAddToRegArith.h"
#include "svmulPackResidualToVRW.h"
#include "svmulSplitToMainAndResidual.h"
#include "unpack.h"
#include "unpackToPartAndCopy.h"
#include "vmmul.h"
#include "vadd.h"
#include "vaddPackResidualToVRW.h"
#include "vaddSplitToMainAndResidual.h"
#include "vrwSVMulAddToRegArith.h"
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

#define DOPARTIALLOOPUNROLLING 1
#define PARTIALUNROLLINGSTARTCOEF 2
#define PARTIALUNROLLINGENDCOEF 32

#if DOCOMPACTLOOPUNROLLING + DOPARTIALLOOPUNROLLING > 1
do you really want to do compact unrolling and partial unrolling?
#endif

void AddGemmTrans() {
  Universe::AddTrans(Gemm::GetClass(), new MMulToMVMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  //Universe::AddTrans(Gemm::GetClass(), new MMulToVMMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(Gemm::GetClass(), new MMulSplitAlongP(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMuSingle, REAL_SINGLE), LLDLALOOPPHASE);

  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMM, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new MMulLoopExp(ABSLAYER, ABSLAYER, DIMK, LLDLAMuDouble, REAL_DOUBLE), LLDLALOOPPHASE);
  
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
  Universe::AddTrans(MVMul::GetClass(), new DistributeMVMulOverMAdd(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuDouble), LLDLALOOPPHASE);
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuDouble), LLDLALOOPPHASE);
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMM, LLDLAMuSingle), LLDLALOOPPHASE);
  Universe::AddTrans(MVMul::GetClass(), new MVMulLoopRef(ABSLAYER, ABSLAYER, DIMN, LLDLAMuSingle), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulToSVMul(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(MVMul::GetClass(), new MVMulToVVDot(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  return;
}

void AddSMMulTrans() {
  Universe::AddTrans(SMMul::GetClass(), new SMulToSVMul(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);
  Universe::AddTrans(SMMul::GetClass(), new SMulToSVMul(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);

  return;
}

void AddUnrollingTrans() {
  if (DOPARTIALLOOPUNROLLING) {
    for (unsigned int mult = PARTIALUNROLLINGSTARTCOEF; mult <= PARTIALUNROLLINGENDCOEF; mult += 2) {
      Universe::AddTrans(SplitSingleIter::GetClass(), new PartiallyUnrollLoop(mult), LLDLALOOPUNROLLPHASE);
    }
  }
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

void AddSVMulAddTrans() {
  Universe::AddTrans(SVMulAdd::GetClass(), new SVMulAddLoopRef(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);
  Universe::AddTrans(SVMulAdd::GetClass(), new VRWSVMulAddToRegArith(ABSLAYER, ABSLAYER), LLDLALOOPPHASE);

  Universe::AddTrans(SVMulAdd::GetClass(), new ResidualSVMulAddToRegArith(ABSLAYER, COLVECTOR), LLDLALOOPPHASE);
  Universe::AddTrans(SVMulAdd::GetClass(), new ResidualSVMulAddToRegArith(ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);

  Universe::AddTrans(SVMulAdd::GetClass(), new SVMulAddSplitToMainAndResidual(ABSLAYER, ABSLAYER, COLVECTOR), LLDLALOOPPHASE);
  Universe::AddTrans(SVMulAdd::GetClass(), new SVMulAddSplitToMainAndResidual(ABSLAYER, ABSLAYER, ROWVECTOR), LLDLALOOPPHASE);
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
  Universe::AddSimp(StoreFromRegs::GetClass(), new EliminateStoreLoadOut(ABSLAYER, ABSLAYER), LLDLARTLPHASE);
  Universe::AddSimp(StoreFromRegs::GetClass(), new EliminateStoreLoad(ABSLAYER, ABSLAYER), LLDLARTLPHASE);
  Universe::AddSimp(StoreFromRegs::GetClass(), new EliminateStore(ABSLAYER, ABSLAYER), LLDLARTLPHASE);

  Universe::AddSimp(UnpackStoreFromRegs::GetClass(), new EliminatePackedStoreLoad(ABSLAYER, ABSLAYER), LLDLARTLPHASE);
  Universe::AddSimp(UnpackStoreFromRegs::GetClass(), new EliminatePackedStore(ABSLAYER, ABSLAYER), LLDLARTLPHASE);

  Universe::AddTrans(Recombine::GetClass(), new EliminateRecombinePartition(ABSLAYER, ABSLAYER), SIMP);
  Universe::AddTrans(Recombine::GetClass(), new EliminateRecombine(ABSLAYER, ABSLAYER), SIMP);

  Universe::AddSimp(LoopTunnel::GetClass(), new HoistDuplicateLoad(), LLDLARTLPHASE);
  Universe::AddSimp(LoopTunnel::GetClass(), new HoistLoadToRegs(), LLDLARTLPHASE);

  Universe::AddTrans(Mul::GetClass(), new AddMulToFMA(ABSLAYER, ABSLAYER), SIMP);
}

void AddPrimPhaseConversions() {
  //  Universe::AddTrans(LoadToRegs::GetClass(), new LoadToContigLoad(ABSLAYER, ABSLAYER), LLDLAPRIMPHASE);
}

void AddArchSpecificTrans() {
  for (auto ext : *arch->SupportedExtensions()) {
    for (auto simp : ext->GetArchTrans()) {
      Universe::AddSimp(simp.first, simp.second, LLDLAPRIMPHASE);
    }
  }
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
  AddSVMulAddTrans();


  AddPrimPhaseConversions();
  AddArchSpecificTrans();
  AddRTLOptimizations();

  AddTransposeTrans();
  AddUnrollingTrans();
  AddSetToZeroTrans();
}

#endif // DOLLDLA
