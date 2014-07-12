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
#include "elemRedist.h"
#include "twoSidedTrxm.h"
#include "blas.h"
#include "loopSupport.h"
#include "helperNodes.h"
#include <time.h>
#include "DLAReg.h"
#include "pack.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lu.h"

#if DOELEM||DOBLIS

#include "driverUtils.h"
#include "debug.h"


//These control which transformations are included
#define USEVR 1
#define USEMCMR 1
#define USESTAR 1
#define USECONTRIB 1
#define USELOCALCOMP 1
#define EXPLOREREDISTS 1
#define EXPLORETRANS 1
#define USESPECIALTRSM 0
#define USELOWERING 1

#define REMOVESCALEBYONE 1 

//good
#define Hetrmm1 1
//bad
#define Hetrmm2 0
#define Hetrmm3 0
//bad
#define TriInv1 0
#define TriInv2 0
//good
#define TriInv3 1
//scalapack
#define TriInv8 0

#if DOELEM
Size smallSize = 500;
Size medSize = 20000;
Size bigSize = 80000;
Size bs = ELEM_BS;
#elif DOSQM || DOSM
Size smallSize = 500;
Size medSize = 8000;
Size bigSize = 10000;
//Size bs = ELEM_BS;
#endif

Trans transA, transB;
Tri tri;
Side side;
char charIn;
Type type;
double size; 
int variant;

PSet* CholExample();
PSet* CholTrsmExample();
PSet* TrsmExample();
PSet* TrmmExample();
PSet* Trmm3Example();
PSet* HegstR1Example();
PSet* HegstR2Example();
PSet* HegstR4Example();
PSet* HegstL1Example();
PSet* HegstL2Example();
PSet* HegstL4Example();
PSet* HegstL5Example();
PSet* HegstLExample();
PSet* HegstRExample();
PSet* GemmExample();
PSet* HemmExample();
PSet* TriInvExample();
PSet* HetrmmExample();
PSet* SPDInvLowerExample();
PSet* PartSPDInvLowerExample();
PSet* SPDInvUpperExample();
PSet* CholHegstExample();
PSet* CholTriInvExample();
//PSet* TriRed();
PSet* Test();
PSet* HerkExample();
PSet* Her2kExample();
PSet* LUExample();
PSet* AppBlkHouseExample();

void AddTrans()
{
  //  MultiTrans *axpyTrans = new MultiTrans;
  //  axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_MC_MR));
  //  Universe::AddTrans(DistAxpy::GetClass(), new DistAxpyToLocalAxpy(D_MC_MR), DPPHASE);
  //  axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_VC_STAR));
  //  axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_STAR_VC));
#if DODPPHASE
  Universe::AddTrans(Axpy::GetClass(), new DistAxpyToLocalAxpy(D_VC_STAR), DPPHASE);
  Universe::AddTrans(Axpy::GetClass(), new DistAxpyToLocalAxpy(D_STAR_VC), DPPHASE);
#if USEVR
  //  axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_VR_STAR));
  //  axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_STAR_VR));
  Universe::AddTrans(Axpy::GetClass(), new DistAxpyToLocalAxpy(D_VR_STAR), DPPHASE);
  Universe::AddTrans(Axpy::GetClass(), new DistAxpyToLocalAxpy(D_STAR_VR), DPPHASE);
#endif
#if USEMCMR
  /*
    axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_MC_STAR));
    axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_STAR_MC));
    axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_MR_STAR));
    axpyTrans->AddTrans(new DistAxpyToLocalAxpy(D_STAR_MR));
  */
  /*
    Universe::AddTrans(Axpy::GetClass(), new DistAxpyToLocalAxpy(D_MC_STAR), DPPHASE);
    Universe::AddTrans(Axpy::GetClass(), new DistAxpyToLocalAxpy(D_STAR_MC), DPPHASE);
    Universe::AddTrans(Axpy::GetClass(), new DistAxpyToLocalAxpy(D_MR_STAR), DPPHASE);
    Universe::AddTrans(Axpy::GetClass(), new DistAxpyToLocalAxpy(D_STAR_MR), DPPHASE);
  */
#endif
  //  Universe::AddTrans(Axpy::GetClass(), axpyTrans, DPPHASE);
#endif //DPPHASE

#if DOSR1PHASE
  Universe::AddTrans(Axpy::GetClass(), new AxpyToBLASAxpy(ABSLAYER, S3LAYER), SR1PHASE);  
#if USELOWERING
  Universe::AddTrans(Axpy::GetClass(), new AxpyLowerLayer(S1LAYER, S3LAYER), SIMP);  
  Universe::AddTrans(Axpy::GetClass(), new AxpyLowerLayer(S2LAYER, S3LAYER), SIMP);  
#endif 
#endif //DOSR1PHASE

#if DODPPHASE
  Universe::AddTrans(Chol::GetClass(), new DistCholToLocalChol, DPPHASE);
  //  Universe::AddTrans(Chol::GetClass(), new CholLoopExp(1), DPPHASE);
  Universe::AddTrans(Chol::GetClass(), new CholLoopExp(2), DPPHASE);
  Universe::AddTrans(Chol::GetClass(), new CholLoopExp(3), DPPHASE);
#endif //DPPHASE

#if DODPPHASE
  Universe::AddTrans(Hetrmm::GetClass(), new DistHetrmmToLocalHetrmm, DPPHASE);
#endif

#if DODPPHASE
#if Hetrmm1
  Universe::AddTrans(Hetrmm::GetClass(), new HetrmmLoopExp(1), DPPHASE);
#endif
#if Hetrmm2
  Universe::AddTrans(Hetrmm::GetClass(), new HetrmmLoopExp(2), DPPHASE);
#endif
#if Hetrmm3
  Universe::AddTrans(Hetrmm::GetClass(), new HetrmmLoopExp(3), DPPHASE);
#endif
#endif //DPPHASE

#if DODPPHASE
  Universe::AddTrans(TriInv::GetClass(), new DistTriInvToLocalTriInv, DPPHASE);

#if TriInv1
  Universe::AddTrans(TriInv::GetClass(), new TriInvLoopExp(1), DPPHASE);
#endif
#if TriInv2
  Universe::AddTrans(TriInv::GetClass(), new TriInvLoopExp(2), DPPHASE);
#endif
#if TriInv3
  Universe::AddTrans(TriInv::GetClass(), new TriInvLoopExp(3), DPPHASE);
#endif
#if TriInv8
  Universe::AddTrans(TriInv::GetClass(), new TriInvLoopExp(8), DPPHASE);
#endif
#endif //DPPHASE

#if DODPPHASE
  Universe::AddTrans(Gemm::GetClass(), new GemmLoopExp(ABSLAYER, DMLAYER, 0), DPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new GemmLoopExp(ABSLAYER, DMLAYER, 1), DPPHASE);
  Universe::AddTrans(Gemm::GetClass(), new GemmLoopExp(ABSLAYER, DMLAYER, 2), DPPHASE);
#endif //DPPHASE

#if DOSR1PHASE
  Universe::AddTrans(Gemm::GetClass(), new GemmLoopExp(ABSLAYER, S1LAYER, 2), SR1PHASE);
#if USELOWERING
  Universe::AddTrans(Gemm::GetClass(), new GemmLowerLayer(ABSLAYER, S1LAYER, DIMN, BLIS_NC_BSVAL), SR1PHASE);
#endif
#endif //DOSR1PHASE
#if DOSR2PHASE
  Universe::AddTrans(Gemm::GetClass(), new SplitGemm(S1LAYER), SR2PHASE);
  Universe::AddTrans(Gemm::GetClass(), new GemmLoopExp(S1LAYER, S2LAYER, 1), SR2PHASE);
  Universe::AddTrans(Gemm::GetClass(), new GemmLoopExp(S1LAYER, S2LAYER, -1), SR2PHASE);
#if USELOWERING
  Universe::AddTrans(Gemm::GetClass(), new GemmLowerLayer(S1LAYER, S2LAYER, DIMK, BLIS_KC_BSVAL), SR2PHASE);
#endif
#endif //SR1PHASE

#if DOSR3PHASE
  Universe::AddTrans(Gemm::GetClass(), new BLISGemmLoopExp(S2LAYER, S3LAYER), SR3PHASE);
#endif //SR2PHASE

#if DODPPHASE  
  MultiTrans *gemmTrans = new MultiTrans;
  gemmTrans->AddTrans(new DistGemmToLocalGemmStatC);
#if USECONTRIB
  gemmTrans->AddTrans(new DistGemmToContribLocalGemmStatANonTrans);
  gemmTrans->AddTrans(new DistGemmToContribLocalGemmStatBNonTrans(false));
  gemmTrans->AddTrans(new DistGemmToContribLocalGemmStatBNonTrans(true));
  gemmTrans->AddTrans(new DistGemmToContribLocalGemmStatBTrans);
  gemmTrans->AddTrans(new DistGemmToContribLocalGemmStatATrans);
#endif
  Universe::AddTrans(Gemm::GetClass(), gemmTrans, DPPHASE);
#endif //DPPHASE

#if DODPPHASE
  Universe::AddTrans(TwoSidedTrxm::GetClass(), new DistTwoSidedTrxmToLocalTwoSidedTrxm, DPPHASE);  
  Universe::AddTrans(TwoSidedTrxm::GetClass(), new TwoSidedTrxmLoopExp(2, ABSLAYER, 
								       TWOSIDEDTRXMCOMPONENTSLAYER, 
								       TWOSIDEDTRXMLAYER), DPPHASE);
  Universe::AddTrans(TwoSidedTrxm::GetClass(), new TwoSidedTrxmLoopExp(4, ABSLAYER, 
								       TWOSIDEDTRXMCOMPONENTSLAYER, 
								       TWOSIDEDTRXMLAYER), DPPHASE);
#endif

#if DOSR1PHASE
  Universe::AddTrans(TwoSidedTrxm::GetClass(), new TwoSidedTrxmLoopExp(2, ABSLAYER, 
								       TWOSIDEDTRXMCOMPONENTSLAYER, 
								       TWOSIDEDTRXMLAYER), SR1PHASE);
  Universe::AddTrans(TwoSidedTrxm::GetClass(), new TwoSidedTrxmLoopExp(4, ABSLAYER, 
								       TWOSIDEDTRXMCOMPONENTSLAYER, 
								       TWOSIDEDTRXMLAYER), SR1PHASE);
#endif

#if DODPPHASE
  Universe::AddTrans(Hemm::GetClass(), new HemmLoopExp(ABSLAYER, DMLAYER, 4), DPPHASE);
  Universe::AddTrans(Hemm::GetClass(), new HemmLoopExp(ABSLAYER, DMLAYER, 8), DPPHASE);
  Universe::AddTrans(Hemm::GetClass(), new HemmLoopExp(ABSLAYER, DMLAYER, -8), DPPHASE);
#endif //DPPHASE


#if DODPPHASE
  MultiTrans *hemmTrans = new MultiTrans;
#if USELOCALCOMP
  hemmTrans->AddTrans(new DistHemmToLocalHemmStatA);
#endif
  hemmTrans->AddTrans(new DistHemmToLocalHemm(D_STAR_VC,D_VC_STAR));
#if USEVR
  hemmTrans->AddTrans(new DistHemmToLocalHemm(D_STAR_VR,D_VR_STAR));
#endif
#if USEMCMR
  hemmTrans->AddTrans(new DistHemmToLocalHemm(D_STAR_MR,D_MR_STAR));
  hemmTrans->AddTrans(new DistHemmToLocalHemm(D_STAR_MC,D_MC_STAR));
#endif
  Universe::AddTrans(Hemm::GetClass(), hemmTrans, DPPHASE);
#endif //DPPHASE

#if DODPPHASE
  MultiTrans *her2kTrans = new MultiTrans;
  her2kTrans->AddTrans(new DistHer2kToLocalTri2k);
#if USECONTRIB
  her2kTrans->AddTrans(new DistHer2kToLocalHer2kContrib(D_STAR_VC,D_VC_STAR));
  her2kTrans->AddTrans(new DistHer2kToLocalHer2kContrib(D_STAR_VR,D_VR_STAR));
#if USEMCMR
  her2kTrans->AddTrans(new DistHer2kToLocalHer2kContrib(D_STAR_MC,D_MC_STAR));
  her2kTrans->AddTrans(new DistHer2kToLocalHer2kContrib(D_STAR_MR,D_MR_STAR));
#endif //USEMCMR
#endif //USECONTRIB
  Universe::AddTrans(Her2k::GetClass(), her2kTrans, DPPHASE);
#endif //DPPHASE

#if DODPPHASE
  Universe::AddTrans(Her2k::GetClass(), new Her2kLoopExp(ABSLAYER, DMLAYER, 1), DPPHASE);
  Universe::AddTrans(Her2k::GetClass(), new Her2kLoopExp(ABSLAYER, DMLAYER, 2), DPPHASE);
  Universe::AddTrans(Her2k::GetClass(), new Her2kLoopExp(ABSLAYER, DMLAYER, 3), DPPHASE);
  Universe::AddTrans(Her2k::GetClass(), new Her2kLoopExp(ABSLAYER, DMLAYER, 4), DPPHASE);
  Universe::AddTrans(Her2k::GetClass(), new Her2kLoopExp(ABSLAYER, DMLAYER, 9), DPPHASE);
#endif

#if DODPPHASE
  Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(ABSLAYER, DMLAYER, 1, LEFT), DPPHASE);
  Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(ABSLAYER, DMLAYER, 2, LEFT), DPPHASE);
  Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(ABSLAYER, DMLAYER, 3, LEFT), DPPHASE);
  Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(ABSLAYER, DMLAYER, 1, RIGHT), DPPHASE);
  Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(ABSLAYER, DMLAYER, 2, RIGHT), DPPHASE);
  Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(ABSLAYER, DMLAYER, 3, RIGHT), DPPHASE);
#endif //DPPHASE

#if DOSR1PHASE
  //  Universe::AddTrans(Trxm::GetClass(), new TrxmRightToLeft<Trxm>(ABSLAYER), DOSR1PHASE);
  Universe::AddTrans(Trmm3::GetClass(), new Trmm3RightToLeft(ABSLAYER), DOSR1PHASE);
#endif //DOSR1PHASE


#if DOSR1PHASE
  Universe::AddTrans(Trxm::GetClass(), new TrmmAxpytoTrxm3(ABSLAYER), SR1PHASE);
  Universe::AddTrans(Trxm::GetClass(), new CopyTrmmtoTrxm3(ABSLAYER), SR1PHASE);
#endif // DOSR1PHASE

#if DOSR1PHASE
  Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(ABSLAYER, S1LAYER, 3, LEFT), SR1PHASE);
  Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(ABSLAYER, S1LAYER, 1, RIGHT), SR1PHASE);
#if USELOWERING
  Universe::AddTrans(Trxm::GetClass(), new TrxmLowerLayer<Trxm>(ABSLAYER, S1LAYER, DIMN, BLIS_NC_BSVAL), SR1PHASE);
#endif

  Universe::AddTrans(Trmm3::GetClass(), new Trmm3LoopExp(ABSLAYER, S1LAYER, 3), SR1PHASE);
#if USELOWERING
    Universe::AddTrans(Trmm3::GetClass(), new TrxmLowerLayer<Trmm3>(ABSLAYER, S1LAYER, DIMN, BLIS_NC_BSVAL), SR1PHASE);
#endif
#endif //SR1PHASE

#if DOSR2PHASE
    Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(S1LAYER, S2LAYER, 2, LEFT), SR2PHASE);
    Universe::AddTrans(Trxm::GetClass(), new TrxmLoopExp(S1LAYER, S2LAYER, 2, RIGHT), SR2PHASE);
#if USELOWERING
  Universe::AddTrans(Trxm::GetClass(), new TrxmLowerLayer<Trxm>(S1LAYER, S2LAYER, DIMK, BLIS_KC_BSVAL), SR2PHASE);
#endif

  Universe::AddTrans(Trmm3::GetClass(), new Trmm3LoopExp(S1LAYER, S2LAYER, 2), SR2PHASE);
#if USELOWERING
  Universe::AddTrans(Trmm3::GetClass(), new TrxmLowerLayer<Trmm3>(S1LAYER, S2LAYER, DIMK, BLIS_KC_BSVAL), SR2PHASE);
#endif
#endif //SR1PHASE

#if DOSR3PHASE
  Universe::AddTrans(Trxm::GetClass(), new BLISTrxmLoopExp(S2LAYER, S3LAYER), SR3PHASE);
  Universe::AddTrans(Trmm3::GetClass(), new BLISTrmm3LoopExp(S2LAYER, S3LAYER), SR3PHASE);
#endif //SR2PHASE


#if DODPPHASE
  Universe::AddTrans(Herk::GetClass(), new DistHerkToLocalTriRK, DPPHASE);
  Universe::AddTrans(Herk::GetClass(), new HerkLoopExp(ABSLAYER, DMLAYER, 1), DPPHASE);
  Universe::AddTrans(Herk::GetClass(), new HerkLoopExp(ABSLAYER, DMLAYER, 2), DPPHASE);
  Universe::AddTrans(Herk::GetClass(), new HerkLoopExp(ABSLAYER, DMLAYER, 5), DPPHASE);
#endif //DODPPHASE



#if DOSR1PHASE
  Universe::AddTrans(TriRK::GetClass(), new TriRKLoopExp(ABSLAYER, S1LAYER, 7), SR1PHASE);
  //  Universe::AddTrans(Herk::GetClass(), new HerkLoopExp(ABSLAYER, SQ1LAYER, 5), SQR1PHASE);
#if USELOWERING
  Universe::AddTrans(TriRK::GetClass(), new TriRKLowerLayer(ABSLAYER, S1LAYER, DIMN, BLIS_NC_BSVAL), SR1PHASE);
#endif
#endif //SQR1PHASE

#if DOSR2PHASE
  Universe::AddTrans(TriRK::GetClass(), new TriRKLoopExp(S1LAYER, S2LAYER, 5), SR2PHASE);
#if USELOWERING
  Universe::AddTrans(TriRK::GetClass(), new TriRKLowerLayer(S1LAYER, S2LAYER, DIMK, BLIS_KC_BSVAL), SR2PHASE);
#endif
#endif //SR2PHASE

#if DOSR3PHASE
  Universe::AddTrans(TriRK::GetClass(), new BLISTriRKLoopExp(S2LAYER, S3LAYER), SR3PHASE);
#endif //SR2PHASE

  
#if DOSR1PHASE
  Universe::AddTrans(Tri2k::GetClass(), new Tri2kLoopExp(ABSLAYER, S1LAYER, 10), SR1PHASE);
#if USELOWERING
  Universe::AddTrans(Tri2k::GetClass(), new Tri2kLowerLayer(ABSLAYER, S1LAYER, DIMN, BLIS_NC_BSVAL), SR1PHASE);
#endif 
#endif //SR1PHASE


#if DOSR2PHASE
  Universe::AddTrans(Tri2k::GetClass(), new Tri2kLoopExp(S1LAYER, S2LAYER, 9), SR1PHASE);
#if USELOWERING
  Universe::AddTrans(Tri2k::GetClass(), new Tri2kLowerLayer(S1LAYER, S2LAYER, DIMK, BLIS_KC_BSVAL), SR2PHASE);
#endif 
#endif //SR2PHASE


#if DOSR3PHASE
  Universe::AddTrans(Tri2k::GetClass(), new Tri2kToTriRK(S2LAYER, S2LAYER), SR3PHASE);
#endif //SQR1PHASE



#if DODPPHASE
  MultiTrans *trxmTrans = new MultiTrans;
  trxmTrans->AddTrans(new DistTrxmToLocalTrxm(D_STAR_VC, D_VC_STAR));
#if USEVR
  trxmTrans->AddTrans(new DistTrxmToLocalTrxm(D_STAR_VR, D_VR_STAR));
#endif
#if USEMCMR
  trxmTrans->AddTrans(new DistTrxmToLocalTrxm(D_STAR_MC, D_MC_STAR));
  trxmTrans->AddTrans(new DistTrxmToLocalTrxm(D_STAR_MR, D_MR_STAR));
#endif
#if USESTAR
  trxmTrans->AddTrans(new DistTrxmToLocalTrxm(D_STAR_STAR, D_STAR_STAR));
#endif
#if USELOCALCOMP
  trxmTrans->AddTrans(new DistTrmmToLocalTrmmStatA);
#endif
  Universe::AddTrans(Trxm::GetClass(), trxmTrans, DPPHASE);

#if USESPECIALTRSM
  Universe::AddTrans(Trxm::GetClass(), new DistTrsmToSpecialLocalTrsm, DPPHASE);
#endif
#endif //DPPHASE



#if DOROPHASE
#if EXPLOREREDISTS
  Universe::AddTrans(RedistNode::GetClass(), new UniqueTransTrans, ROPHASE);


  Universe::AddTrans(RedistNode::GetClass(), new ExpandRedistribution<D_MC_MR,D_MC_STAR>, ROPHASE);
  Universe::AddTrans(RedistNode::GetClass(), new ExpandRedistribution<D_MC_MR,D_STAR_MR>, ROPHASE);
  Universe::AddTrans(RedistNode::GetClass(), new ExpandRedistribution<D_MC_MR,D_MR_STAR>, ROPHASE);
  Universe::AddTrans(RedistNode::GetClass(), new ExpandRedistribution<D_MC_MR,D_STAR_MC>, ROPHASE);
  Universe::AddTrans(RedistNode::GetClass(), new ExpandRedistribution<D_VC_STAR,D_MR_STAR>, ROPHASE);

#endif
#endif 

#if DOROPHASE
#if EXPLORETRANS

  for (int i = 0; i < 4; ++i) {
    Universe::AddTrans(Tri2k::GetClass(), new Tri2kTrans(i, TRANS), ROPHASE);
    Universe::AddTrans(Tri2k::GetClass(), new Tri2kTrans(i, CONJTRANS), ROPHASE);
  }

  for (int i = 0; i < 2; ++i) {
    Universe::AddTrans(Gemm::GetClass(), new GemmTrans(i, TRANS), ROPHASE);
    Universe::AddTrans(Gemm::GetClass(), new GemmTrans(i, CONJTRANS), ROPHASE);
    Universe::AddTrans(TriRK::GetClass(), new TriRKTrans(i, TRANS), ROPHASE);
    Universe::AddTrans(TriRK::GetClass(), new TriRKTrans(i, CONJTRANS), ROPHASE);    
  }

  Universe::AddTrans(RedistNode::GetClass(), new RedistTrans(TRANS), ROPHASE);    
  Universe::AddTrans(RedistNode::GetClass(), new RedistTrans(CONJTRANS), ROPHASE);    

  //Universe::AddTrans(LocalGemv::GetClass(), new GemvTrans(0, TRANS), ROPHASE);    
  
  Universe::AddTrans(Trxm::GetClass(), new TrxmTrans(TRANS), ROPHASE);
  Universe::AddTrans(Trxm::GetClass(), new TrxmTrans(CONJTRANS), ROPHASE);
#endif
#endif

#if DOELEM
  GemmInputReordering *t0 = new GemmInputReordering(0);
  GemmInputReordering *t1 = new GemmInputReordering(1);
  t0->m_inverse = t1;
  t1->m_inverse = t0;
#endif

#if DODPPHASE
  Universe::AddTrans(Gemm::GetClass(), t0, DPPHASE);
  Universe::AddTrans(Gemm::GetClass(), t1, DPPHASE);
#endif //DODPPHASE

#if DOSR1PHASE
  Universe::AddTrans(LU::GetClass(), new LULoopExp(ABSLAYER, ABSLAYER, S3LAYER, 5), SR1PHASE);
#endif //DOSR1PHASE

#if DOSMPPHASE
  //Universe::AddTrans(Split::GetClass(), new IncreaseParallelizedLoop, SMPPHASE);

#if NUMPROCS>1
  Universe::AddTrans(LoopTunnel::GetClass(), new ParallelizeOuterNDim(ALLPROCCOMM), SMPPHASE);
    //Universe::AddTrans(LoopTunnel::GetClass(), new ParallelizeOuterNDim(ALLPROCCOMM), SIMP);
#endif //NUMPROCS>1

#if NUMPROCS>1
  //  Universe::AddTrans(LoopTunnel::GetClass(), new ParallelizeK(ALLPROCCOMM), SMPPHASE);
  //  Universe::AddTrans(Axpy::GetClass(), new ParallelizeAxpy(S3LAYER, ALLPROCCOMM), SMPPHASE);
#endif //NUMPROCS>1

#if NUML2PERPROC>1
  Universe::AddTrans(PackBuff::GetClass(), new ParallelizeMDim(PROCCOMM), SMPPHASE);
  Universe::AddTrans(PackBuff::GetClass(), new ParallelizeMDim(ALLL2COMM), SMPPHASE);
#endif //NUML2PERPROC>1

#if NUMCORESPERL2>1
  Universe::AddTrans(PackBuff::GetClass(), new ParallelizeInnerNDim(PROCCOMM), SMPPHASE);
  Universe::AddTrans(PackBuff::GetClass(), new ParallelizeInnerNDim(L2COMM), SMPPHASE);
  Universe::AddTrans(PackBuff::GetClass(), new ParallelizeInnerNDim(L2COMMSUBALLL2), SMPPHASE);
  Universe::AddTrans(PackBuff::GetClass(), new ParallelizeInnerNDim(ALLL2COMM), SMPPHASE);
#endif //NUMCORESPERL2>1

#endif


}

void AddSimplifiers()
{ 
#if REMOVESCALEBYONE
  Universe::AddTrans(ScaleNode::GetClass(), new RemoveScaleByOne, SIMP);
#endif //REMOVESCALEBYONE

#if DOELEM
  Universe::AddTrans(RedistNode::GetClass(), 
		     new UseTransposedRedist(D_MC_MR,D_STAR_MC_T,D_STAR_MC_H), SIMP);

  Universe::AddTrans(RedistNode::GetClass(), 
		     new UseTransposedRedist(D_MR_MC,D_STAR_MR_T,D_STAR_MR_H), SIMP);
 
  Universe::AddTrans(RedistNode::GetClass(), new ExpandRedistribution<D_MC_MR,D_VR_STAR>, SIMP);
  Universe::AddTrans(Trxm::GetClass(), new LTrmmToTrsm, SIMP);
  Universe::AddTrans(Trxm::GetClass(), new DTrmmToTrsm, SIMP);
  for(int i=1; i < D_LASTDIST; ++i) {
    Universe::AddTrans(RedistNode::GetClass(), new RemoveNOPRedistribs((DistType)i), SIMP);
    Universe::AddTrans(RedistNode::GetClass(), new RemoveWastedRedist((DistType)i), SIMP);
    for(int j=1; j < D_LASTDIST; ++j)
      if(i!=j) {
        Universe::AddTrans(RedistNode::GetClass(), new CombineRedistribs((DistType)i,(DistType)j), SIMP);
      }
  }
  //Universe::AddTrans(PackB::GetClass(), new CombinePackB, SIMP);
  Universe::AddTrans(MakeTrapNode::GetClass(), new MoveMakeTrap, SIMP);

  //The four wings
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_VC_STAR, D_MC_STAR, D_MC_MR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_MC_STAR, D_VC_STAR), SIMP);

  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_STAR_VR, D_STAR_MR, D_MC_MR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_MR, D_STAR_VR), SIMP);

  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_VR_STAR, D_MR_STAR, D_MR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MR_MC, D_MR_STAR, D_VR_STAR), SIMP);

  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_STAR_VC, D_STAR_MC, D_MR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MR_MC, D_STAR_MC, D_STAR_VC), SIMP);

  //A transpose op
  //  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_VC_STAR, D_STAR_MC_T, D_MC_MR), SIMP);

  //Through the upper middle
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_VC_STAR, D_MR_STAR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MR_STAR, D_VC_STAR, D_MC_MR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_VR_STAR, D_MR_STAR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MR_STAR, D_VR_STAR, D_MC_MR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_VC_STAR, D_MR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MR_MC, D_VC_STAR, D_MC_MR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_VR_STAR, D_MR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MR_MC, D_VR_STAR, D_MC_MR), SIMP);

  /*
  //Through the bottom middle
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_VR, D_STAR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_STAR_MC, D_STAR_VR, D_MC_MR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_VC, D_STAR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_STAR_MC, D_STAR_VC, D_MC_MR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_VR, D_MR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MR_MC, D_STAR_VR, D_MC_MR), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_VC, D_MR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MR_MC, D_STAR_VC, D_MC_MR), SIMP);
  */

  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_VR, D_STAR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_VC, D_STAR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_VR, D_MR_MC), SIMP);
  Universe::AddTrans(RedistNode::GetClass(), new FindMidDistributions(D_MC_MR, D_STAR_VC, D_MR_MC), SIMP);

#endif

#if DOSR1PHASE
  Universe::AddTrans(Hemm::GetClass(), new BLISHemmToGemm(ABSLAYER), SIMP);
#endif

#if DOSR2PHASE
  Universe::AddTrans(PackBuff::GetClass(), new LoopInvariantPackBuffMotion, GLOBSIMP);
  Universe::AddTrans(PackBuff::GetClass(), new UnifyPackBuffParams, SIMP);
  Universe::AddTrans(Pack::GetClass(), new LoopInvariantPackMotion, GLOBSIMP);
  Universe::AddTrans(Pack::GetClass(), new CombinePacking, GLOBSIMP);
  Universe::AddTrans(Pack::GetClass(), new ReuseTrsmPacking(S3LAYER), SIMP);
  Universe::AddTrans(PackBuff::GetClass(), new CombinePackBuff, SIMP);
  Universe::AddTrans(Transpose::GetClass(), new CombineTranspose, SIMP);
#endif //SR2PHASE

#if DOSR1PHASE
  Universe::AddTrans(Herk::GetClass(), new HerkToTriRK(ABSLAYER), SIMP);
  Universe::AddTrans(Her2k::GetClass(), new Her2kToTri2K(ABSLAYER), SIMP);
#endif //DOSR1PHASE

#if DOSOPHASE
  Universe::AddTrans(PackBuff::GetClass(), new RenamePackBuff, SIMP);
#endif //DOSOPHASE

}

void Usage()
{
  cout << "./driver arg1 arg2 arg3 arg4\n";
  cout <<" arg1 == 0  -> Load from file arg1\n";
  cout <<"         1  -> Chol L/U\n";
  cout <<"         2  -> Chol + TriInv L\n";
  cout <<"         3  -> Trsm 0/1/2 L/R L/U N/T\n";
  cout <<"         4  -> Hegst Lower var/(0 for all variants) L/R\n";
  cout <<"         5  -> Gemm A/B/C N/T N/T\n";
  cout <<"         6  -> Tri Inv + Hetrmm\n";
  cout <<"         7  -> Chol + Hegst side\n";
  cout <<"         8  -> Special Test\n";
  cout <<"         9  -> Triangular matrix inverse L/U\n";
  cout <<"         10 -> Hetrmm L/U\n";
  cout <<"         11 -> SPD inversion (lower)\n";
  cout <<"         12 -> SPD inversion (upper)\n";
  cout <<"         13 -> Chol - both variants\n";
  cout <<"         14 -> Chol + Trsm L\n";
  cout <<"         15 -> Herk C/R L/U N/T/C\n";
  cout <<"         16 -> Her2k C/R L/U N/T/C\n";
  cout <<"         17 -> Hemm A/C C/R L/R L/U\n";
  cout <<"         18 -> Trmm A/C L/R L/U N/T\n";
  cout <<"         19 -> LU\n";
  cout <<"         20 -> Apply Block Householder\n";
  cout <<"         21 -> Trmm3 A/C L/R L/U N/T\n";
  cout <<" arg2 == 0 -> PrintAll\n";
  //  cout <<"         0  -> Print optimal graph number\n";
  cout <<"        num -> Print graph number num\n";
  cout <<" arg3 == -1 -> iterate until completion\n";
  cout <<"        num -> take num iterations\n";
}

int main(int argc, const char* argv[])
{
  //    omp_set_num_threads(1);
#ifdef _OPENMP
  omp_set_nested(true);
#endif
  //  PrintType printType = CODE;
  int numIters = -1;
  PSet* (*algFunc)();
  //  unsigned int whichGraph = 0;
  int algNum;
  string fileName;

  if(argc < 2) {
    Usage();
    return 0;
  }
  else {
    algNum = atoi(argv[1]);
    switch(algNum) {
    case (0):
      fileName = argv[2];
    case (1):
      algFunc = CholExample;
      tri = CharToTri(*argv[2]);
      break;
    case (2):
      algFunc = CholTriInvExample;
      break;
    case (3):
      algFunc = TrsmExample;
      size = atoi(argv[2]);
      side = CharToSide(*argv[3]);
      tri = CharToTri(*argv[4]);
      transA = CharToTrans(*argv[5]);
      break;
    case (4):
      variant = atoi(argv[2]);
      side = CharToSide(*argv[3]);
      if (side == RIGHT) {
	if (variant == 0) {
	  algFunc = HegstRExample;
	}
	else if (variant == 1) {
	  algFunc = HegstR1Example;
	}
	else if (variant == 2) {
	  algFunc = HegstR2Example;
	}
	else if (variant == 4) {
	  algFunc = HegstR4Example;
	}
	else
	  throw;
      }
      else
	switch (variant){
	case(0):
	  algFunc = HegstLExample;
	  break;
	case(1):
	  algFunc = HegstL1Example;
	  break;
	case(2):
	  algFunc = HegstL2Example;
	  break;
	case(4):
	  algFunc = HegstL4Example;
	  break;
	case(5):
	  algFunc = HegstL5Example;
	  break;
	default:
	  throw;
	}
      break;
    case (5):
      algFunc = GemmExample;
      charIn = *argv[2];
      transA = CharToTrans(*argv[3]);
      transB = CharToTrans(*argv[4]);
      break;
    case (6):
      algFunc = PartSPDInvLowerExample;
      break;
    case (7):
      algFunc = CholHegstExample;
      side = CharToSide(*argv[2]);
      break;
    case (8):
      algFunc = Test;
      break;
    case (9):
      algFunc = TriInvExample;
      tri = CharToTri(*argv[2]);
      break;
    case (10):
      algFunc = HetrmmExample;
      tri = CharToTri(*argv[2]);
      break;
    case (11):
      algFunc = SPDInvLowerExample;
      break;
    case (12):
      algFunc = SPDInvUpperExample;
      break;
    case (13):
      throw;
      break;
    case (14):
      algFunc = CholTrsmExample;
      break;
    case (15):
      algFunc = HerkExample;
      type = CharToType(*argv[2]);
      tri = CharToTri(*(argv[3]));
      transA = CharToTrans(*(argv[4]));
      break;
    case (16):
      algFunc = Her2kExample;
      type = CharToType(*argv[2]);
      tri = CharToTri(*(argv[3]));
      transA = CharToTrans(*(argv[4]));
      break;
    case (17):
      algFunc = HemmExample;
      charIn = *(argv[2]);
      type = CharToType(*argv[3]);
      side = CharToSide(*argv[4]);
      tri = CharToTri(*(argv[5]));
      break;
    case (18):
      algFunc = TrmmExample;
      charIn = *(argv[2]);
      side = CharToSide(*argv[3]);
      tri = CharToTri(*argv[4]);
      transA = CharToTrans(*argv[5]);
      break;
    case (19):
      algFunc = LUExample;
      break;
    case (20):
      algFunc = AppBlkHouseExample;
      break;
    case (21):
      algFunc = Trmm3Example;
      charIn = *(argv[2]);
      side = CharToSide(*argv[3]);
      tri = CharToTri(*argv[4]);
      transA = CharToTrans(*argv[5]);
      break;
    default:
      Usage();
      return 0;
    }
  }

  RegAllDLANodes();
  AddTrans();
  AddSimplifiers();

  Universe uni;
  time_t start, start2, end;
  uni.PrintStats();

  if (algNum==0) {
    time(&start);
    uni.Init(fileName);
    time(&end);
    cout << "Unflatten took " << difftime(end,start) << " seconds\n";
    //    uni.SanityCheck();
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
  else {
    uni.Init(algFunc());
    time(&start);
  }


#if DODPPHASE
  if (CurrPhase == DPPHASE) {
    cout << "Expanding DP phase\n";
    uni.Expand(-1, DPPHASE, DLACullDP);
    time(&end);
    cout << "DP phase took " << difftime(end,start) << " seconds\n";

    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOROPHASE
  if (CurrPhase == ROPHASE) {
    cout << "Expanding RO phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, ROPHASE, DLACullRO);
    time(&end);
    cout << "RO phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOSR1PHASE
  if (CurrPhase == SR1PHASE) {
    cout << "Expanding SR1 phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SR1PHASE, DLACullSR);
    time(&end);
    cout << "SR1 phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOSR2PHASE
  if (CurrPhase == SR2PHASE) {
    cout << "Expanding SR2 phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SR2PHASE, DLACullSR);
    time(&end);
    cout << "SR2 phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif


#if DOSR3PHASE
  if (CurrPhase == SR3PHASE) {
    cout << "Expanding SR3 phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SR3PHASE, DLACullSR);
    time(&end);
    cout << "SR3 phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif


#if DOSOPHASE
  if (CurrPhase == SOPHASE) {
    cout << "Expanding SO phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SOPHASE, DLACullSR);
    time(&end);
    cout << "SO phase took " << difftime(end,start2) << " seconds\n";
    
    cout << "Propagating\n";
    cout.flush();
    time(&start2);
    uni.Prop();
    time(&end);
    cout << "Propagation took " << difftime(end,start2) << " seconds\n";
  }
#endif

#if DOSMPPHASE
  if (CurrPhase == SMPPHASE) {
    cout << "Shared-memory parallelization phase\n";
    cout << "Starting with " << uni.TotalCount() << endl;
    time(&start2);
    uni.Expand(numIters, SMPPHASE, DLACullSR);
    time(&end);
    cout << "SMP phase took " << difftime(end,start2) << " seconds\n";
    
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

  uni.PrintAll(algNum);

  /*  if (whichGraph <= 0)
      uni.PrintAll();
      else
      uni.Print(cout, CODE, whichGraph); */

  return 0;
}


PSet* TrsmExample()
{
  InputNode *Lin;
  InputNode *Xin;

  if (side == LEFT) {
    if (size == 0) {
      Lin = new InputNode("L input", smallSize, smallSize, tri==LOWER?"L":"U");
      Xin = new InputNode("X input", smallSize, smallSize, "X");
    }
    else if (size == 1) {
      Lin = new InputNode("L input", smallSize, smallSize, tri==LOWER?"L":"U");
      Xin = new InputNode("X input", smallSize, smallSize, "X");
    }
    else if (size == 2) {
      Lin = new InputNode("L input", medSize, medSize, tri==LOWER?"L":"U");
      Xin = new InputNode("X input", medSize, bigSize, "X");
    }
    else
      throw;
  }
  else {
    if (size == 0) {
      Lin = new InputNode("L input", smallSize, smallSize, tri==LOWER?"L":"U");
      Xin = new InputNode("X input", smallSize, smallSize, "X");
    }
    else if (size == 1) {
      Lin = new InputNode("L input", smallSize, smallSize, tri==LOWER?"L":"U");
      Xin = new InputNode("X input", smallSize, smallSize, "X");
    }
    else if (size == 2) {
      Lin = new InputNode("L input", medSize, medSize, tri==LOWER?"L":"U");
      Xin = new InputNode("X input", bigSize, medSize, "X");
    }
    else
      throw;
  }

  Trxm *loop = new Trxm(true, ABSLAYER, side, tri, NONUNIT, transA, COEFONE, REAL);
  loop->AddInputs(4,
		  Lin, 0,
		  Xin, 0);
			       
  OutputNode *Bout = new OutputNode("B output");
  Bout->AddInput(loop,0);

  Poss *poss2 = new Poss(Bout,true);
  PSet *set2 = new PSet (poss2);

  return set2;

}

PSet* CholExample()
{
#if DOELEM
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  PossTunnel *tun = new PossTunnel(POSSTUNIN);
  tun->AddInput(Ain,0);

  Chol *loop = new Chol(ABSLAYER, tri);
  loop->AddInput(tun,0);

  Poss *innerPoss = new Poss(loop,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
#else
  throw;
#endif
}

PSet* CholTrsmExample()
{
#if DOELEM
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");
  InputNode *Bin = new InputNode("B input", bigSize, bigSize, "B");

  PossTunnel *tun = new PossTunnel(POSSTUNIN);
  tun->AddInput(Ain,0);

  PossTunnel *tun2 = new PossTunnel(POSSTUNIN);
  tun2->AddInput(Bin,0);

  Chol *loop = new Chol(ABSLAYER, tri);
  loop->AddInput(tun,0);

  Trxm *trsm = new Trxm(true, ABSLAYER, LEFT, LOWER, NONUNIT, NORMAL, COEFONE, COMPLEX);
  trsm->AddInput(loop);
  trsm->AddInput(tun2);

  Poss *innerPoss = new Poss(trsm,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
#else
  throw;
#endif
  
}

PSet* CholTriInvExample()
{
#if DOELEM
  InputNode *Ain = new InputNode("L input", bigSize, bigSize, "L");

  PossTunnel *tun = new PossTunnel(POSSTUNIN);
  tun->AddInput(Ain,0);

  Chol *loop = new Chol(ABSLAYER, LOWER);
  loop->AddInput(tun,0);

  TriInv *loop2 = new TriInv(ABSLAYER, LOWER);
  loop2->AddInput(loop,0);

  Poss *innerPoss = new Poss(loop2,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
#else
  throw;
#endif
}


PSet* HegstRExample()
{
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

  TwoSidedTrxm *hegst = new TwoSidedTrxm(ABSLAYER, true, LOWER);
  hegst->AddInputs(4,
		   Lin, 0,
		   Ain, 0);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(hegst, 0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}


PSet* HegstR1Example()
{
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

  Loop *loop = TwoSidedTrsmLowerVar1Alg(Lin, 0,
				      Ain, 0,
				      TWOSIDEDTRXMCOMPONENTSLAYER, TWOSIDEDTRXMLAYER);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(loop->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout, true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}


PSet* HegstR4Example()
{
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

  Loop *loop = TwoSidedTrsmLowerVar4Alg(Lin, 0,
				      Ain, 0,
				      TWOSIDEDTRXMCOMPONENTSLAYER, TWOSIDEDTRXMLAYER);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(loop->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout, true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}

PSet* HegstR2Example()
{
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

  Loop *loop = TwoSidedTrsmLowerVar2Alg(Lin, 0,
					Ain, 0,
					TWOSIDEDTRXMCOMPONENTSLAYER, TWOSIDEDTRXMLAYER);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(loop->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout, true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}


PSet* HegstL1Example()
{
#if DOELEM
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain,0);
  splitA->SetUpStats(FULLUP, NOTUP,
		     FULLUP, NOTUP);

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");
  SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitL->AddInput(Lin,0);
  splitL->SetAllStats(FULLUP);

#if DOELEM
  TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y21");
  throw; //SetLayer?
#else
  TempVarNode *Yin = new TempVarNode("Y21");
  Yin->SetLayer(S3LAYER);
#endif
  Yin->AddInput(splitA, 5);

  Hemm *hemm = new Hemm(DMLAYER, LEFT, LOWER, COEFONE, COEFZERO, COMPLEX);
  hemm->AddInput(splitA,8);
  hemm->AddInput(splitL,5);
  hemm->AddInput(Yin,0);
  PSet *set1 = new PSet(new Poss(hemm,false));

  Trxm *trmm = new Trxm(false, DMLAYER, LEFT, LOWER, NONUNIT, NORMAL, COEFONE, COMPLEX);
  trmm->AddInputs(4, splitL, 4, splitA, 5);
  Poss *poss2 = new Poss(trmm,false);
  PSet *set2 = new PSet(poss2);
  Axpy *axpy1 = new Axpy(DMLAYER, COEFONEHALF);
  axpy1->AddInput(set1->OutTun(0),0);
  axpy1->AddInput(set2->OutTun(0),0);
  PSet *set6 = new PSet(new Poss(axpy1,false));


  TwoSidedTrxm *hegst = new TwoSidedTrxm(DMLAYER, false, LOWER);
  hegst->AddInput(splitL,4);
  hegst->AddInput(splitA,4);

  PSet *set4 = new PSet(new Poss(hegst,false));


  Her2k *her2k = new Her2k(DMLAYER, LOWER, CONJTRANS, COEFONE, COEFONE, COMPLEX);
  her2k->AddInput(set6->OutTun(0),0);
  her2k->AddInput(splitL,5);
  her2k->AddInput(set4->OutTun(0),0);
  PSet *set7 = new PSet(new Poss(her2k,false));

  Axpy *axpy2 = new Axpy(DMLAYER, COEFONEHALF);
  axpy2->AddInput(set1->OutTun(0),0);
  axpy2->AddInput(set6->OutTun(0),0);
  PSet *set8 = new PSet(new Poss(axpy2,false));

  Trxm *trmm2 = new Trxm(false, DMLAYER, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm2->AddInput(splitL,8);
  trmm2->AddInput(set8->OutTun(0),0);
  PSet *set5 = new PSet(new Poss(trmm2,false));

  CombineSingleIter *comA;
  comA = splitA->CreateMatchingCombine(2, 
				       4, set7->OutTun(0), 0,
				       5, set5->OutTun(0), 0);

  CombineSingleIter *comL;
  comL = splitL->CreateMatchingCombine(0);

  Poss *loopPoss = new Poss(2,
			    comA,
			    comL);

  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(loop->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;

#else
  throw;
#endif
}


PSet* HegstL2Example()
{
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

  Loop *loop = TwoSidedTrmmLowerVar2Alg(Lin, 0,
					    Ain, 0,
					    TWOSIDEDTRXMCOMPONENTSLAYER, TWOSIDEDTRXMLAYER);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(loop->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}


PSet* HegstL4Example()
{
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

  Loop *loop = TwoSidedTrmmLowerVar4Alg(Lin, 0, Ain, 0,
				     TWOSIDEDTRXMCOMPONENTSLAYER, TWOSIDEDTRXMLAYER);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(loop->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}

PSet* HegstLExample()
{
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

  TwoSidedTrxm *hegst = new TwoSidedTrxm(ABSLAYER, false, LOWER);
  hegst->AddInputs(4, 
		   Lin, 0,
		   Ain, 0);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(hegst, 0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}


PSet* HegstL5Example()
{
#if DOELEM
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");
  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain,0);
  splitA->SetUpStats(PARTUP, FULLUP,
		     PARTUP, NOTUP);

  InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");
  SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitL->AddInput(Lin,0);
  splitL->SetAllStats(FULLUP);

#if DOELEM
  TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y10");
  throw; //setlayer ?
#else
  TempVarNode *Yin = new TempVarNode("Y10");
  Yin->SetLayer(S3LAYER);
#endif
  Yin->AddInput(splitA, 1);

  Hemm *hemm = new Hemm(DMLAYER, LEFT, LOWER, COEFONE, COEFZERO, COMPLEX);
  hemm->AddInput(splitA,4);
  hemm->AddInput(splitL,1);
  hemm->AddInput(Yin,0);
  PSet *set1 = new PSet(new Poss(hemm,false));

  Trxm *trmm = new Trxm(false, DMLAYER, RIGHT, LOWER, NONUNIT, NORMAL, COEFONE, COMPLEX);
  trmm->AddInputs(4, splitL, 0, splitA, 1);
  Poss *poss2 = new Poss(trmm,false);
  PSet *set2 = new PSet(poss2);

  Axpy *axpy1 = new Axpy(DMLAYER, COEFONEHALF);
  axpy1->AddInput(set1->OutTun(0),0);
  axpy1->AddInput(set2->OutTun(0),0);
  PSet *set6 = new PSet(new Poss(axpy1,false));

  Her2k *her2k = new Her2k(DMLAYER, LOWER, CONJTRANS, COEFONE, COEFONE, COMPLEX);
  her2k->AddInput(splitL,1);
  her2k->AddInput(set6->OutTun(0),0);
  her2k->AddInput(splitA,0);
  PSet *set7 = new PSet(new Poss(her2k,false));

  Axpy *axpy2 = new Axpy(DMLAYER, COEFONEHALF);
  axpy2->AddInput(set1->OutTun(0),0);
  axpy2->AddInput(set6->OutTun(0),0);
  PSet *set8 = new PSet(new Poss(axpy2,false));
  
  Trxm *trmm2 = new Trxm(false, DMLAYER, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm2->AddInput(splitL,4);
  trmm2->AddInput(set8->OutTun(0),0);
  PSet *set5 = new PSet(new Poss(trmm2,false));

  TwoSidedTrxm *hegst = new TwoSidedTrxm(DMLAYER, false, LOWER);
  hegst->AddInput(splitL,4);
  hegst->AddInput(splitA,4);
  PSet *set4 = new PSet(new Poss(hegst,false));

  CombineSingleIter *comA;
  comA = splitA->CreateMatchingCombine(3,
				       0, set7->OutTun(0), 0,
				       1, set5->OutTun(0), 0,
				       4, set4->OutTun(0), 0);
				       

  CombineSingleIter *comL;
  comL = splitL->CreateMatchingCombine(0);
  
  Poss *loopPoss = new Poss(2,
			    comA,
			    comL);
#if DOELEM
  Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(loop->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
#else
  throw;
#endif
#endif
}

PSet* GemmExample()
{
  InputNode *Ain;
  InputNode *Bin;
  InputNode *Cin;

  switch(charIn)
    {
    case('A'):
      Ain = new InputNode("A input", bigSize, bigSize, "A");
      if (transB == NORMAL)
	Bin = new InputNode("B input", bigSize, medSize, "B");
      else
	Bin = new InputNode("B input", medSize, bigSize, "B");
      Cin = new InputNode("C input", bigSize, medSize, "C");
      break;
    case('B'):
      if (transA == NORMAL)
	Ain = new InputNode("A input", medSize, bigSize, "A");
      else
	Ain = new InputNode("A input", bigSize, medSize, "A");
      Bin = new InputNode("B input", bigSize, bigSize, "B");
      Cin = new InputNode("C input", medSize, bigSize, "C");
      break;     
    case('C'):
      if (transA == NORMAL)
	Ain = new InputNode("A input", bigSize, medSize, "A");
      else
	Ain = new InputNode("A input", medSize, bigSize, "A");
      if (transB == NORMAL)
	Bin = new InputNode("B input", medSize, bigSize, "B");
      else
	Bin = new InputNode("B input", bigSize, medSize, "B");
      Cin = new InputNode("C input", bigSize, bigSize, "C");
      break;
    default:
      throw;
    }

  Gemm *loop = new Gemm(ABSLAYER, transA, transB, COEFONE, COEFONE, REAL);
  loop->AddInputs(6,
		  Ain, 0,
		  Bin, 0,
		  Cin, 0);
			       
  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(loop,0);

  Poss *poss2 = new Poss(Cout,true);
  PSet *set2 = new PSet (poss2);

  return set2;
}



PSet* TriInvExample()
{
  InputNode *Ain = new InputNode("input", bigSize, bigSize, tri==LOWER ? "L" : "U");

  PossTunnel *tun = new PossTunnel(POSSTUNIN);
  tun->AddInput(Ain,0);

  TriInv *loop = new TriInv(ABSLAYER, tri);
  loop->AddInput(tun,0);

  Poss *innerPoss = new Poss(loop,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}

PSet* HetrmmExample()
{
#if DOELEM
  InputNode *Ain = new InputNode("input", bigSize, bigSize, tri == LOWER ? "L" : "U");

  PossTunnel *tun = new PossTunnel(POSSTUNIN);
  tun->AddInput(Ain,0);

  Hetrmm *loop = new Hetrmm(ABSLAYER, tri);
  loop->AddInput(tun,0);

  Poss *innerPoss = new Poss(loop,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
#else
  throw;
#endif
}

PSet* PartSPDInvLowerExample()
{
#if DOELEM
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  PossTunnel *tun = new PossTunnel(POSSTUNIN);
  tun->AddInput(Ain,0);

  TriInv *loop2 = new TriInv(ABSLAYER, LOWER);
  loop2->AddInput(tun,0);

  Hetrmm *loop3 = new Hetrmm(ABSLAYER, LOWER);
  loop3->AddInput(loop2,0);

  Poss *innerPoss = new Poss(loop3,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
#else
  throw;
#endif
}

PSet* SPDInvLowerExample()
{
#if DOELEM
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  PossTunnel *tun = new PossTunnel(POSSTUNIN);
  tun->AddInput(Ain,0);

  Chol *loop = new Chol(ABSLAYER, LOWER);
  loop->AddInput(tun,0);

  TriInv *loop2 = new TriInv(ABSLAYER, LOWER);
  loop2->AddInput(loop,0);

  Hetrmm *loop3 = new Hetrmm(ABSLAYER, LOWER);
  loop3->AddInput(loop2,0);

  Poss *innerPoss = new Poss(loop3,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
#else
  throw;
#endif
}

PSet* SPDInvUpperExample()
{
#if DOELEM
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");

  PossTunnel *tun = new PossTunnel(POSSTUNIN);
  tun->AddInput(Ain,0);

  Chol *loop = new Chol(ABSLAYER, UPPER);
  loop->AddInput(tun,0);

  TriInv *loop2 = new TriInv(ABSLAYER, UPPER);
  loop2->AddInput(loop,0);

  Hetrmm *loop3 = new Hetrmm(ABSLAYER, UPPER);
  loop3->AddInput(loop2,0);
  
  Poss *innerPoss = new Poss(loop3,true);
  PSet *innerSet = new PSet(innerPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(innerSet->OutTun(0),0);

  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
#else
  throw;
#endif
}


PSet* Test()
{
  return NULL;
}


/*
  PSet* TriRed()
  {
  InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");
  InputNode *Yin = new InputNode("Y input", bigSize, bigSize, "Y");
  InputNode *Uin = new InputNode("U input", bigSize, bigSize, "U");

  SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitA->AddInput(Ain,0);
  splitA->SetUpStats(FULLUP, FULLUP,
  FULLUP, PARTUP);

  SplitSingleIter *splitU = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitU->AddInput(Uin,0);
  splitU->SetUpStats(FULLUP, FULLUP,
  FULLUP, PARTUP);

  SplitSingleIter *splitY = new SplitSingleIter(PARTDIAG, POSSTUNIN);
  splitY->AddInput(Yin,0);
  splitY->SetUpStats(FULLUP, FULLUP,
  FULLUP, PARTUP);

  //alpha11 := alpha11 - 2 * u10^T * y10
  
  InputNode *temp1 = new InputNode("temp1", ONE, ONE, "temp1");

  DistDot *dot1 = new DistDot(-2);
  dot1->AddInputs(6, 
  splitU, 1,
  splitY, 1,
  temp1, 0);
  Poss *poss1 = new Poss(dot1, false);
  PSet *set1 = new PSet(poss1);

  ConstVal *const1 = new ConstVal("constVal", COEFONE);

  Axpy *axpy1 = new Axpy(DMLAYER);
  axpy1->AddInputs(6,
  const1, 0,
  set1->OutTun(0), 0,
  splitA, 4);
  Poss *poss2 = new Poss(axpy1, false);
  PSet *set2 = new PSet(poss2);

  //a21 := a21 - U20 * y10 - Y20 * u10

  DistGemv *gemv1 = new DistGemv(NORMAL, -1, 1);
  gemv1->AddInputs(6,
  splitU, 2,
  splitY, 1,
  splitA, 5);
  Poss *poss3 = new Poss(gemv1, false);
  PSet *set3 = new PSet(poss3);

  DistGemv *gemv2 = new DistGemv(NORMAL, -1, 1);
  gemv2->AddInputs(6,
  splitY, 2,
  splitU, 1,
  set3->OutTun(0), 0);
  Poss *poss4 = new Poss(gemv2, false);
  PSet *set4 = new PSet(poss4);

  //[u21, tau, a21] := HOUSEV( a21 )

  InputNode *tau = new InputNode("tau", ONE, ONE, "tau");
  DistHouseV *house = new DistHouseV;
  house->AddInputs(6, 
  splitU, 5,
  tau, 0,
  set4->OutTun(0), 0);
  Poss *poss5 = new Poss(house,false);
  PSet *set5 = new PSet(poss5);

  //y21 := A22 * u21
  DistGemv *gemv3 = new DistGemv(NORMAL, 1, 0);
  gemv3->AddInputs(6,
  splitA, 8,
  set5->OutTun(0), 0,
  splitY, 5);
  Poss *poss6 = new Poss(gemv3,false);
  PSet *set6 = new PSet(poss6);

  //temp2 := U20^T * u21
  InputNode *temp2 = new InputNode("temp2", ((DLANode*)set6->OutTun(0))->GetM(0), ONE, "temp2");
  DistGemv *gemv4 = new DistGemv(TRANS, 1, 0);
  gemv4->AddInputs(6,
  splitU, 2,
  set5->OutTun(0), 0,
  temp2, 0);
  Poss *poss7 = new Poss(gemv4, false);
  PSet *set7 = new PSet(poss7);

  //y21 := y21 - Y20 * (U20^T * u21)
  DistGemv *gemv5 = new DistGemv(NORMAL, -1, 1);
  gemv5->AddInputs(6,
  splitY, 2,
  set7->OutTun(0), 0,
  set6->OutTun(0), 0);
  Poss *poss8 = new Poss(gemv5, false);
  PSet *set8 = new PSet(poss8);

  //temp3 := Y20^T * u21
  InputNode *temp3 = new InputNode("temp3", ((DLANode*)set8->OutTun(0))->GetM(0), ONE, "temp3");
  Poss *tempPoss = new Poss(temp3, false);
  PSet *tempSet = new PSet(tempPoss);

  DistGemv *gemv6 = new DistGemv(TRANS, 1, 0);
  gemv6->AddInputs(6,
  splitY, 2,
  set5->OutTun(0), 0,
  tempSet->OutTun(0), 0);
  Poss *poss9 = new Poss(gemv6, false);
  PSet *set9 = new PSet(poss9);

  //y21 := y21 - Y20 * (U20^T * u21) - U20 * (Y20^T * u21)
  DistGemv *gemv7 = new DistGemv(NORMAL, -1, 1);
  gemv7->AddInputs(6,
  splitU, 2,
  set9->OutTun(0), 0,
  set8->OutTun(0), 0);
  Poss *poss10 = new Poss(gemv7, false);
  PSet *set10 = new PSet(poss10);
  

  // beta := -1 * ut2^T * y21 / 2
  InputNode *beta = new InputNode("beta", ONE, ONE, "beta");
  tempPoss = new Poss(beta, false);
  tempSet = new PSet(tempPoss);

  DistDot *dot2 = new DistDot(COEFNEGONEHALF);
  dot2->AddInputs(6, 
  set5->OutTun(0), 0,
  set10->OutTun(0), 0,
  tempSet->OutTun(0), 0);
  Poss *poss11 = new Poss(dot2, false);
  PSet *set11 = new PSet(poss11);

  // beta := beta / tau
  ScalInvert *inver1 = new ScalInvert;
  inver1->AddInput( set5->OutTun(1), 0);
  Poss *invertPoss = new Poss(inver1, false);
  PSet *invertSet = new PSet(invertPoss);

  DistScal *scal1 = new DistScal;
  scal1->AddInputs(4,
  invertSet->OutTun(0), 0,
  set11->OutTun(0), 0);
  Poss *poss12 = new Poss(scal1, false);
  PSet *set12 = new PSet(poss12);

  RedistNode *redist = new RedistNode(D_STAR_STAR);
  redist->AddInput(set12->OutTun(0), 0);
  tempPoss = new Poss(redist);
  tempSet = new PSet(tempPoss);

  //y21 := y21 - beta * u21 / tau
  Axpy *axpy2 = new Axpy(DMLAYER);
  axpy2->AddInputs(6,
  tempSet->OutTun(0), 0,
  set5->OutTun(0), 0,
  set10->OutTun(0), 0);
  Poss *poss13 = new Poss(axpy2, false);
  PSet *set13 = new PSet(poss13);

  //y21 := y21 / tau
  DistScal *scal2 = new DistScal;
  scal2->AddInputs(4,
  invertSet->OutTun(0), 0,
  set13->OutTun(0), 0);
  Poss *poss14 = new Poss(scal2, false);
  PSet *set14 = new PSet(poss14);		   

  CombineSingleIter *comA = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comA->AddInput(splitA,0);
  comA->AddInput(splitA,1);
  comA->AddInput(splitA,2);
  comA->AddInput(splitA,3);
  comA->AddInput(set2->OutTun(0),0);
  comA->AddInput(set5->OutTun(2),0);
  comA->AddInput(splitA,6);
  comA->AddInput(splitA,7);
  //comA->AddInput(set->OutTun(0),0);
  //  comA->AddInput(set1->OutTun(0),0);
  comA->AddInput(splitA,8);
  //comA->AddInput(set2->OutTun(0),0);
  comA->AddInput(splitA,9);

  comA->CopyUpStats(splitA);

  CombineSingleIter *comU = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comU->AddInput(splitU,0);
  comU->AddInput(splitU,1);
  comU->AddInput(splitU,2);
  comU->AddInput(splitU,3);
  comU->AddInput(splitU,4);
  //  comU->AddInput(splitU,5);
  comU->AddInput(set5->OutTun(0),0);
  comU->AddInput(splitU,6);
  comU->AddInput(splitU,7);
  comU->AddInput(splitU,8);
  comU->AddInput(splitU,9);

  comU->CopyUpStats(splitU);

  CombineSingleIter *comY = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
  comY->AddInput(splitY,0);
  comY->AddInput(splitY,1);
  comY->AddInput(splitY,2);
  comY->AddInput(splitY,3);
  comY->AddInput(splitY,4);
  //  comY->AddInput(splitY,5);
  comY->AddInput(set14->OutTun(0),0);
  comY->AddInput(splitY,6);
  comY->AddInput(splitY,7);
  comY->AddInput(splitY,8);
  comY->AddInput(splitY,9);

  comY->CopyUpStats(splitY);

  Poss *loopPoss = new Poss(3, comA, comU, comY);
  Loop *loop = new Loop(ELEMLOOP, loopPoss);

  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(loop->OutTun(0),0);

  OutputNode *Uout = new OutputNode("U output");
  Uout->AddInput(loop->OutTun(1),0);

  OutputNode *Yout = new OutputNode("Y output");
  Yout->AddInput(loop->OutTun(2),0);
   
  Poss *outerPoss = new Poss(Aout,true);
  PSet *outerSet = new PSet(outerPoss);

  return outerSet;
  }*/

PSet* HerkExample()
{
  InputNode *Ain;
#if DODM
  if (transA == NORMAL)
    Ain = new InputNode("A input", bigSize, medSize, "A", D_MC_MR);
  else
    Ain = new InputNode("A input", medSize, bigSize, "A", D_MC_MR);

  InputNode *Cin = new InputNode("C input", bigSize, bigSize, "C", D_MC_MR);
#else
  if (transA == NORMAL)
    Ain = new InputNode("A input", bigSize, medSize, "A");
  else
    Ain = new InputNode("A input", medSize, bigSize, "A");

  InputNode *Cin = new InputNode("C input", bigSize, bigSize, "C");
#endif

  Herk *loop = new Herk(ABSLAYER, tri, transA, COEFONE, COEFONE, type);
  loop->AddInputs(4,
		  Ain, 0,
		  Cin, 0);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(loop,0);

  Poss *poss = new Poss(Cout,true);
  PSet *set = new PSet(poss);
  
  return set;
}

PSet* Her2kExample()
{
  InputNode *Ain, *Bin;
#if DOELEM
  if (transA == NORMAL) {
    Ain = new InputNode("A input", bigSize, medSize, "A", D_MC_MR);
    Bin = new InputNode("B input", bigSize, medSize, "B", D_MC_MR);
  }
  else {
    Ain = new InputNode("A input", medSize, bigSize, "A", D_MC_MR);
    Bin = new InputNode("B input", medSize, bigSize, "B", D_MC_MR);
  }

  InputNode *Cin = new InputNode("C input", bigSize, bigSize, "C", D_MC_MR);
#else
  if (transA == NORMAL) {
    Ain = new InputNode("A input", bigSize, medSize, "A");
    Bin = new InputNode("B input", bigSize, medSize, "B");
  }
  else {
    Ain = new InputNode("A input", medSize, bigSize, "A");
    Bin = new InputNode("B input", medSize, bigSize, "B");
  }

  InputNode *Cin = new InputNode("C input", bigSize, bigSize, "C");
#endif

  Her2k *loop = new Her2k(ABSLAYER, tri, transA, COEFONE, COEFONE, type);
  loop->AddInputs(6,
		  Ain, 0,
		  Bin, 0,
		  Cin, 0);

  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(loop,0);

  Poss *poss = new Poss(Cout,true);
  PSet *set = new PSet(poss);
  
  return set;
}


PSet* HemmExample()
{
  InputNode *Ain;
  InputNode *Bin;
  InputNode *Cin;

  switch(charIn)
    {
    case('A'):
      Ain = new InputNode("A input", medSize, medSize, "A");
      if (side == LEFT) {
	Bin = new InputNode("B input", medSize, smallSize, "B");
	Cin = new InputNode("C input", medSize, smallSize, "C");
      }
      else {
	Bin = new InputNode("B input", smallSize, medSize, "B");
	Cin = new InputNode("C input", smallSize, medSize, "C");
      }
      break;
    case('C'):
      Ain = new InputNode("A input", medSize, medSize, "A");
      if (side == LEFT) {
	Bin = new InputNode("B input", medSize, bigSize, "B");	
	Cin = new InputNode("C input", medSize, bigSize, "C");
      }
      else {
	Bin = new InputNode("B input", bigSize, medSize, "B");
	Cin = new InputNode("C input", bigSize, medSize, "C");
      }
      break;
    default:
      throw;
    }

  Hemm *loop = new Hemm(ABSLAYER, side, tri, COEFONE, COEFONE, type);
  loop->AddInputs(6,
		  Ain, 0,
		  Bin, 0,
		  Cin, 0);
			       
  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(loop,0);

  Poss *poss2 = new Poss(Cout,true);
  PSet *set2 = new PSet (poss2);

  return set2;
}



PSet* TrmmExample()
{
  InputNode *Lin;
  InputNode *Xin;

  if (side == LEFT) {
    switch(charIn)
      {
      case('A'):
	Lin = new InputNode("L input", bigSize, bigSize, tri==LOWER?"L":"U");
	Xin = new InputNode("X input", bigSize, medSize, "X");
	break;
      case('C'):
	Lin = new InputNode("L input", bigSize, bigSize, tri==LOWER?"L":"U");
	Xin = new InputNode("X input", bigSize, bigSize, "X");
	break;
      default:
	throw;
      }
  }
  else {
    switch(charIn)
      {
      case('A'):
	Lin = new InputNode("L input", bigSize, bigSize, tri==LOWER?"L":"U");
	Xin = new InputNode("X input", medSize, bigSize, "X");
	break;
      case('C'):
	Lin = new InputNode("L input", bigSize, bigSize, tri==LOWER?"L":"U");
	Xin = new InputNode("X input", bigSize, bigSize, "X");
	break;
      default:
	throw;
      }
  }

  Trxm *loop = new Trxm(false, ABSLAYER, side, tri, NONUNIT, transA, COEFONE, REAL);
  loop->AddInputs(4,
		  Lin, 0,
		  Xin, 0);
			       
  OutputNode *Bout = new OutputNode("B output");
  Bout->AddInput(loop,0);

  Poss *poss2 = new Poss(Bout,true);
  PSet *set2 = new PSet (poss2);

  return set2;

}

PSet* Trmm3Example()
{
  InputNode *Lin;
  InputNode *Xin;
  InputNode *Cin;
  

  if (side == LEFT) {
    switch(charIn)
      {
      case('A'):
	Lin = new InputNode("L input", medSize, medSize, tri==LOWER?"L":"U");
	Xin = new InputNode("X input", medSize, smallSize, "X");
	Cin = new InputNode("C input", medSize, smallSize, "C");
	break;
      case('C'):
	Lin = new InputNode("L input", medSize, medSize, tri==LOWER?"L":"U");
	Xin = new InputNode("X input", medSize, bigSize, "X");
	Cin = new InputNode("C input", medSize, bigSize, "C");
	break;
      default:
	throw;
      }
  }
  else {
    switch(charIn)
      {
      case('A'):
	Lin = new InputNode("L input", medSize, medSize, tri==LOWER?"L":"U");
	Xin = new InputNode("X input", smallSize, medSize, "X");
	Cin = new InputNode("C input", medSize, medSize, "C");
	break;
      case('C'):
	Lin = new InputNode("L input", medSize, medSize, tri==LOWER?"L":"U");
	Xin = new InputNode("X input", bigSize, medSize, "X");
	Cin = new InputNode("C input", bigSize, medSize, "C");
	break;
      default:
	throw;
      }
  }

  Trmm3 *loop = new Trmm3(ABSLAYER, side, tri, NONUNIT, transA, COEFONE, COEFONE, REAL);
  loop->AddInputs(6,
		  Lin, 0,
		  Xin, 0,
		  Cin, 0);
			       
  OutputNode *Cout = new OutputNode("C output");
  Cout->AddInput(loop,0);

  Poss *poss2 = new Poss(Cout,true);
  PSet *set2 = new PSet (poss2);

  return set2;
}

PSet* CholHegstExample()
{
#if DOELEM
  if (side == RIGHT) {
    InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");
    SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
    splitA->AddInput(Ain,0);
    splitA->SetUpStats(PARTUP, FULLUP,
		       PARTUP, PARTUP);

    InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

    Chol *cholloop = new Chol(ABSLAYER, LOWER);
    cholloop->AddInput(Lin,0);

    SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
    splitL->AddInput(cholloop,0);
    splitL->SetAllStats(FULLUP);

    TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y21");
    Yin->SetLayer(ABSLAYER);
    Yin->AddInput(splitA, 5);

    Trxm *trsm = new Trxm(true, DMLAYER, LEFT, LOWER, NONUNIT, NORMAL, COEFONE,COMPLEX);
    trsm->AddInput(splitL,4);
    trsm->AddInput(splitA,1);
    PSet *set1 = new PSet(new Poss(trsm,false));

    TwoSidedTrxm *hegst = new TwoSidedTrxm(DMLAYER, true, LOWER);
    hegst->AddInput(splitL,4);
    hegst->AddInput(splitA,4);
    PSet *set2 = new PSet(new Poss(hegst,false));

    Gemm *gemm = new Gemm(DMLAYER, NORMAL, NORMAL, COEFNEGONE, COEFONE, COMPLEX);
    gemm->AddInput(splitL,5);
    gemm->AddInput(set1->OutTun(0),0);
    gemm->AddInput(splitA,2);
    PSet *set3 = new PSet(new Poss(gemm,false));

    Hemm *hemm = new Hemm(DMLAYER, RIGHT, LOWER, COEFONEHALF, COEFZERO, COMPLEX);
    hemm->AddInput(set2->OutTun(0),0);
    hemm->AddInput(splitL,5);
    hemm->AddInput(Yin,0);
    PSet *set4 = new PSet(new Poss(hemm,false));

    Trxm *Trsm2 = new Trxm(true, DMLAYER, RIGHT, LOWER, NONUNIT, CONJTRANS,COEFONE,COMPLEX);
    Trsm2->AddInput(splitL,4);
    Trsm2->AddInput(splitA,5);
    PSet *set5 = new PSet(new Poss(Trsm2,false));


    Axpy *axpy1 = new Axpy(DMLAYER, COEFONEHALF);
    axpy1->AddInput(set4->OutTun(0),0);
    axpy1->AddInput(set5->OutTun(0),0);
    PSet *set6 = new PSet(new Poss(axpy1,false));
  
    Her2k *her2k = new Her2k(DMLAYER, LOWER, NORMAL, COEFNEGONE, COEFONE, COMPLEX);
    her2k->AddInput(splitL,5);
    her2k->AddInput(set6->OutTun(0),0);
    her2k->AddInput(splitA,8);
    PSet *set7 = new PSet(new Poss(her2k,false));

    Axpy *axpy2 = new Axpy(DMLAYER, COEFONEHALF);
    axpy2->AddInput(set4->OutTun(0),0);
    axpy2->AddInput(set6->OutTun(0),0);
    PSet *set8 = new PSet(new Poss(axpy2,false));

    CombineSingleIter *comA = new CombineSingleIter(PARTDIAG, POSSTUNOUT);
    comA->AddInput(splitA,0);
    comA->AddInput(set1->OutTun(0),0);
    comA->AddInput(set3->OutTun(0),0);
    comA->AddInput(splitA,3);
    comA->AddInput(set2->OutTun(0),0);
    comA->AddInput(set8->OutTun(0),0);
    comA->AddInput(splitA,6);
    comA->AddInput(splitA,7);
    comA->AddInput(set7->OutTun(0),0);
    comA->AddInput(splitA,9);
  
    comA->CopyTunnelInfo(splitA);

    CombineSingleIter *comL = splitL->CreateMatchingCombine(0);

    Poss *loopPoss = new Poss(2,
			      comA,
			      comL);
#if DOELEM
    Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

    OutputNode *Aout = new OutputNode("A output");
    Aout->AddInput(loop->OutTun(0),0);

    OutputNode *Lout = new OutputNode("L output");
    Lout->AddInput(cholloop,0);


    Poss *outerPoss = new Poss(2, Aout, Lout);
    PSet *outerSet = new PSet(outerPoss);
  
    return outerSet;
#else
    throw;
#endif
  }
  else {
    InputNode *Ain = new InputNode("A input", bigSize, bigSize, "A");
    SplitSingleIter *splitA = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
    splitA->AddInput(Ain,0);
    splitA->SetUpStats(PARTUP, FULLUP,
		       PARTUP, NOTUP);

    InputNode *Lin = new InputNode("L input", bigSize, bigSize, "L");

    Chol *cholloop = new Chol(ABSLAYER, LOWER);
    cholloop->AddInput(Lin,0);

    SplitSingleIter *splitL = new SplitSingleIter(PARTDIAG, POSSTUNIN);
    splitL->AddInput(cholloop,0);
    splitL->SetAllStats(FULLUP);

#if DOELEM
    TempVarNode *Yin = new TempVarNode(D_MC_MR, "Y10");
    Yin->SetLayer(SMLAYER);
#else
    TempVarNode *Yin = new TempVarNode("Y10");
    Yin->SetLayer(S3LAYER);
#endif
    Yin->AddInput(splitA, 1);

    Hemm *hemm = new Hemm(DMLAYER, LEFT, LOWER, COEFONE, COEFZERO, COMPLEX);
    hemm->AddInput(splitA,4);
    hemm->AddInput(splitL,1);
    hemm->AddInput(Yin,0);
    PSet *set1 = new PSet(new Poss(hemm,false));

    Axpy *axpy1 = new Axpy(DMLAYER, COEFONEHALF);
    axpy1->AddInput(set1->OutTun(0),0);
    axpy1->AddInput(splitA,1);
    PSet *set6 = new PSet(new Poss(axpy1,false));

    Her2k *her2k = new Her2k(DMLAYER, LOWER, CONJTRANS, COEFONE, COEFONE, COMPLEX);
    her2k->AddInput(set6->OutTun(0),0);
    her2k->AddInput(splitL,1);
    her2k->AddInput(splitA,0);
    PSet *set7 = new PSet(new Poss(her2k,false));

    Axpy *axpy2 = new Axpy(DMLAYER, COEFONEHALF);
    axpy2->AddInput(set1->OutTun(0),0);
    axpy2->AddInput(set6->OutTun(0),0);
    PSet *set8 = new PSet(new Poss(axpy2,false));

    Trxm *trmm = new Trxm(false, DMLAYER, LEFT, LOWER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
    trmm->AddInputs(4, splitL, 4, set8->OutTun(0), 0);
    Poss *poss2 = new Poss(trmm,false);
    PSet *set2 = new PSet(poss2);
  
    TwoSidedTrxm *hegst = new TwoSidedTrxm(DMLAYER, false, LOWER);
    hegst->AddInput(splitL,4);
    hegst->AddInput(splitA,4);
    PSet *set4 = new PSet(new Poss(hegst,false));

    Gemm *gemm = new Gemm(DMLAYER, NORMAL, NORMAL, COEFONE, COEFONE, COMPLEX);
    gemm->AddInput(splitA, 5);
    gemm->AddInput(splitL, 1);
    gemm->AddInput(splitA, 2);
    PSet *set3 = new PSet(new Poss(gemm,false));

    Trxm *trmm2 = new Trxm(false, DMLAYER, RIGHT, LOWER, NONUNIT, NORMAL, COEFONE, COMPLEX);
    trmm2->AddInput(splitL,4);
    trmm2->AddInput(splitA,5);
    PSet *set5 = new PSet(new Poss(trmm2,false));

    CombineSingleIter *comA;
    comA = splitA->CreateMatchingCombine(5,
					 0,set7->OutTun(0),0,
					 1,set2->OutTun(0),0,
					 2,set3->OutTun(0),0,
					 4,set4->OutTun(0),0,
					 5,set5->OutTun(0),0);

    CombineSingleIter *comL = splitL->CreateMatchingCombine(0);

    Poss *loopPoss = new Poss(2,
			      comA,
			      comL);
#if DOELEM
    Loop *loop = new Loop(ELEMLOOP, loopPoss, USEELEMBS);

    OutputNode *Aout = new OutputNode("A output");
    Aout->AddInput(loop->OutTun(0),0);

    OutputNode *Lout = new OutputNode("L output");
    Lout->AddInput(cholloop,0);

    Poss *outerPoss = new Poss(2, Aout, Lout);
    PSet *outerSet = new PSet(outerPoss);
  
    return outerSet;
#else
    throw;
#endif
  }
#else
  throw;
#endif
}

PSet* LUExample()
{
  InputNode *A = new InputNode("A input", bigSize, bigSize, "A");
  InputNode *P = new InputNode("P input", bigSize, 1, "P");

  LU *lu = new LU(ABSLAYER);
  lu->AddInputs(4,
		A, 0,
		P, 0);
			       
  OutputNode *Aout = new OutputNode("A output");
  Aout->AddInput(lu,0);

  OutputNode *Pout = new OutputNode("P output");
  Pout->AddInput(lu,1);

  Poss *poss2 = new Poss(2, Aout, Pout);
  PSet *set2 = new PSet (poss2);

  return set2;
}

PSet* AppBlkHouseExample()
#if DOBLIS
{
  InputNode *U = new InputNode("U input", BLIS_KC_BSVAL*10, BLIS_KC_BSVAL*10, "U");
  InputNode *T = new InputNode("T input", BLIS_KC_BSVAL, BLIS_KC_BSVAL*10, "T");
  InputNode *B = new InputNode("B input", BLIS_KC_BSVAL*10, BLIS_KC_BSVAL*10, "B");
  
  SplitSingleIter *splitU = new SplitSingleIter(PARTDIAG, POSSTUNIN, true);
  splitU->AddInput(U,0);
  splitU->SetAllStats(FULLUP);

  SplitSingleIter *splitB = new SplitSingleIter(PARTDOWN, POSSTUNIN);
  splitB->AddInput(B,0);
  splitB->SetUpStats(FULLUP, FULLUP,
		     PARTUP, PARTUP);

  SplitSingleIter *splitT = new SplitSingleIter(PARTRIGHT, POSSTUNIN);
  splitT->AddInput(T,0);
  splitT->SetAllStats(FULLUP);

  ViewTL *view = new ViewTL(S3LAYER);
  view->AddInputs(6,
		  splitT, 1,
		  splitT, 1,
		  splitB, 1);


  TempVarNode *W = new TempVarNode("W");

  W->SetLayer(S3LAYER);
  W->AddInput(splitB, 1);

  ViewTL *view2 = new ViewTL(S3LAYER);
  view2->AddInputs(6,
		   W, 0,
		   splitB, 1,
		   splitB, 1);

  Copy *copy = new Copy;
  copy->SetLayer(ABSLAYER);
  copy->AddInputs(4, 
		  splitB, 1,
		  view2, 0);

  Trxm *trmm1 = new Trxm(false, ABSLAYER, LEFT, LOWER, UNIT, CONJTRANS, COEFONE, COMPLEX);
  trmm1->AddInputs(4,
		   splitU, 4,
		   copy, 0);

  Gemm *gemm1 = new Gemm(ABSLAYER, CONJTRANS, NORMAL, COEFONE, COEFONE, COMPLEX);
  gemm1->AddInputs(6,
		   splitU, 5,
		   splitB, 2,
		   trmm1, 0);
  
  Trxm *trsm1 = new Trxm(true, ABSLAYER, LEFT, UPPER, NONUNIT, CONJTRANS, COEFONE, COMPLEX);
  trsm1->AddInputs(4,
		   view, 0,
		   gemm1, 0);

  Gemm *gemm2 = new Gemm(ABSLAYER, NORMAL, NORMAL, COEFNEGONE, COEFONE, COMPLEX);
  gemm2->AddInputs(6,
		   splitU, 5,
		   trsm1, 0,
		   splitB, 2);

  Trxm *trmm2 = new Trxm(false, ABSLAYER, LEFT, LOWER, UNIT, NORMAL, COEFNEGONE, COMPLEX);
  trmm2->AddInputs(4,
		   splitU, 4,
		   trsm1, 0);

  Axpy *axpy1 = new Axpy(ABSLAYER, COEFONE);
  axpy1->AddInputs(4,
		   trmm2, 0,
		   splitB, 1);

  CombineSingleIter *comB = splitB->CreateMatchingCombine(2,
						1, axpy1, 0,
						2, gemm2, 0);
  
  CombineSingleIter *comU = splitU->CreateMatchingCombine(0);

  CombineSingleIter *comT = splitT->CreateMatchingCombine(0);
  
  Poss *loopPoss = new Poss(3, comB, comU, comT);
  Loop *loop = new Loop(BLISLOOP, loopPoss, USEBLISOUTERBS);

  OutputNode *Bout = new OutputNode("B output");
  Bout->AddInput(loop->OutTun(0),0);

  Poss *outerPoss = new Poss(1, Bout);
  PSet *outerSet = new PSet(outerPoss);
  
  return outerSet;
}
#else
{
  throw;
}
#endif

#endif //DOELEM || DOBLIS
