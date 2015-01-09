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

#include "ccsd.h"

#if DOTENSORS

#include "helperNodes.h"
#include "tensorHelpers.h"
#include "yaxppx.h"
#include "zaxpby.h"
#include "yaxpy.h"
#include "yxpy.h"
#include "zaxpbyppx.h"
#include "contraction.h"

InputNode *CreateInput2(string name, Size size1, Size size2)
{
    Sizes sizes[2];
    sizes[0].AddRepeatedSizes(size1,1,1);
    sizes[1].AddRepeatedSizes(size2,1,1);
    return new InputNode(name, sizes, name, 2);
}

InputNode *CreateInput4(string name, Size size1, Size size2, Size size3, Size size4)
{
    Sizes sizes[4];
    sizes[0].AddRepeatedSizes(size1,1,1);
    sizes[1].AddRepeatedSizes(size2,1,1);
    sizes[2].AddRepeatedSizes(size3,1,1);
    sizes[3].AddRepeatedSizes(size4,1,1);
    return new InputNode(name, sizes, name, 4);
}

RealPSet* W_bmje_calc(DLANode *w_bmje, DLANode *x_bmej,
		      DLANode *r_bmef, DLANode *t_fj,
		      DLANode *u_mnje, DLANode *v_femn,
		      DLANode *T_bfnj,
		      DLANode *Tau_efim,
		      const Size big, const Size small)
{
  InputNode *W_bmje = CreateInput4("W_bmje", big, small, small, big);


  TempVarNode *temp1 = new TempVarNode(T_bfnj->DataType(0).GetDist(), "temp1");
  temp1->AddInput(T_bfnj,0);
  TempVarNode *temp2 = new TempVarNode(v_femn->DataType(0).GetDist(), "temp2");
  temp2->AddInput(v_femn,0);
  TempVarNode *temp3 = new TempVarNode(u_mnje->DataType(0).GetDist(), "temp3");
  temp3->AddInput(u_mnje,0);
  TempVarNode *temp4 = new TempVarNode(r_bmef->DataType(0).GetDist(), "temp4");
  temp4->AddInput(r_bmef,0);

/*  
  InputNode *temp1 = CreateInput4("temp1", big, big, small, small);
  InputNode *temp2 = CreateInput4("temp2", big, big, small, small);
  InputNode *temp3 = CreateInput4("temp3", small, small, small, big);
  InputNode *temp4 = CreateInput4("temp4", big, small, big, big);
*/

  YAxpPx *axpy0 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONEHALF, "bmej", "bmje");
  axpy0->AddInputs0(3,
		     w_bmje, x_bmej, W_bmje);
  Poss *axpy0Poss = new Poss(axpy0);
  RealPSet * axpy0Set = new RealPSet(axpy0Poss);
  
  ZAxpBypPx *axpy1 = new ZAxpBypPx(ABSLAYER, COEFONEHALF, COEFNEGONE, "bfjn", "bfnj");
  axpy1->AddInputs0(4,
		    T_bfnj, Tau_efim, T_bfnj, temp1);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);
  

  YAxpPx *axpy2 = new YAxpPx(ABSLAYER, COEFNEGONE, COEFTWO, "femn", "fenm");
  axpy2->AddInputs0(3,
		     v_femn, v_femn, temp2);
  Poss *axpy2Poss = new Poss(axpy2);
  RealPSet * axpy2Set = new RealPSet(axpy2Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"femn","bfnj","bmje",(string)"fn");
  cont2->AddInputs0(3,
		    axpy2Set->OutTun(0), axpy1Set->OutTun(0), axpy0Set->OutTun(0));
  Poss *cont2Poss = new Poss(cont2);
  RealPSet *cont2Set = new RealPSet(cont2Poss);


  YAxpPx *axpy3 = new YAxpPx(ABSLAYER, COEFNEGONE, COEFTWO, "mnje", "nmje");
  axpy3->AddInputs0(3,
		    u_mnje, u_mnje, temp3);
  Poss *axpy3Poss = new Poss(axpy3);
  RealPSet * axpy3Set = new RealPSet(axpy3Poss);

  Contraction *cont3 = new Contraction(ABSLAYER,COEFNEGONE,COEFONE,REAL,"mnje","bn","bmje",(string)"n");
  cont3->AddInputs0(3,
		    axpy3Set->OutTun(0), t_fj, cont2Set->OutTun(0));
  Poss *cont3Poss = new Poss(cont3);
  RealPSet *cont3Set = new RealPSet(cont3Poss);


  YAxpPx *axpy4 = new YAxpPx(ABSLAYER, COEFNEGONE, COEFTWO, "bmef", "bmfe");
  axpy4->AddInputs0(3,
		    r_bmef, r_bmef, temp4);
  Poss *axpy4Poss = new Poss(axpy4);
  RealPSet * axpy4Set = new RealPSet(axpy4Poss);

  Contraction *cont4 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"bmef","fj","bmje",(string)"f");
  cont4->AddInputs0(3,
		    axpy4Set->OutTun(0), t_fj, cont3Set->OutTun(0));
  Poss *cont4Poss = new Poss(cont4);
  RealPSet *cont4Set = new RealPSet(cont4Poss);

  return cont4Set;
}

RealPSet* X_bmej_calc(DLANode *x_bmej, DLANode *r_bmef,
		      DLANode *t_fj, 
		      DLANode *u_mnje, DLANode *v_femn,
		      DLANode *T_bfnj,
		      DLANode *Tau_efim,
		      const Size big, const Size small)
{
  InputNode *X_bmej = CreateInput4("X_bmej", big, small, big, small);


  TempVarNode *temp1 = new TempVarNode(T_bfnj->DataType(0).GetDist(), "temp1");
  temp1->AddInput(T_bfnj,0);
  //  InputNode *temp1 = CreateInput4("temp1", big, big, small, small);

  Copy *copy1 = new Copy(SMLAYER);
  copy1->AddInputs0(2,
		   x_bmej,
		   X_bmej);
  Poss *copy1Poss = new Poss(copy1);
  RealPSet *copy1Set = new RealPSet(copy1Poss);

  ZAxpBy *axpy1 = new ZAxpBy(SMLAYER, COEFONE, COEFNEGONE);
  axpy1->AddInputs0(3,
		    Tau_efim, T_bfnj,
		    temp1);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFNEGONE,COEFZERO,REAL,"femn","bfnj","bmej",(string)"fn");
  cont2->AddInputs0(3,
		    v_femn, axpy1Set->OutTun(0), copy1Set->OutTun(0));
  Poss *cont2Poss = new Poss(cont2);
  RealPSet *cont2Set = new RealPSet(cont2Poss);


  Contraction *cont3 = new Contraction(ABSLAYER,COEFNEGONE,COEFONE,REAL,"mnje","bn","bmej",(string)"n");
  cont3->AddInputs0(3,
		    u_mnje, t_fj, cont2Set->OutTun(0));
  Poss *cont3Poss = new Poss(cont3);
  RealPSet *cont3Set = new RealPSet(cont3Poss);

  Contraction *cont4 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"bmef","fj","bmej",(string)"f");
  cont4->AddInputs0(3,
		    r_bmef, t_fj, cont3Set->OutTun(0));
  Poss *cont4Poss = new Poss(cont4);
  RealPSet *cont4Set = new RealPSet(cont4Poss);

  return cont4Set;
}

RealPSet* U_mnie_calc(DLANode *t_fj, 
		      DLANode *u_mnie, DLANode *v_femn,
		      const Size big, const Size small)
{
  InputNode *U_mnie = CreateInput4("U_mnie", small, small, small, big);


  Copy *copy1 = new Copy(SMLAYER);
  copy1->AddInputs0(2,
		   u_mnie,
		   U_mnie);
  Poss *copy1Poss = new Poss(copy1);
  RealPSet *copy1Set = new RealPSet(copy1Poss);

  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"femn","fi","mnie",(string)"f");
  cont1->AddInputs0(3,
		    v_femn, t_fj, copy1Set->OutTun(0));
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  return cont1Set;
}

RealPSet* Q_mnij_calc(DLANode *q_mnij,
		      DLANode *t_fj, 
		      DLANode *u_mnie, DLANode *v_efmn,
		      DLANode *T_efij,
		      DLANode *Tau_efim,
		      const Size big, const Size small)
{
  InputNode *Q_mnij = CreateInput4("Q_mnij", small, small, small, small);

  TempVarNode *temp1 = new TempVarNode(Q_mnij->DataType(0).GetDist(), "temp1");
  temp1->AddInput(Q_mnij,0);
  /*
  InputNode *temp1 = CreateInput4("temp1", big, big, small, small);
  InputNode *temp2 = CreateInput4("temp2", small, small, small, small);
  */

  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"mnie","ej","mnij",(string)"e");
  cont1->AddInputs0(3,
		    u_mnie, t_fj, temp1);
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  
  YAxpPx *axpy1 = new YAxpPx(ABSLAYER, COEFONE, COEFONE, "mnij", "nmji");
  axpy1->AddInputs0(3,
		    cont1Set->OutTun(0), cont1Set->OutTun(0), Q_mnij);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);

  Contraction *cont3 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"efmn","efij","mnij",(string)"ef");
  cont3->AddInputs0(3,
		    v_efmn, Tau_efim, axpy1Set->OutTun(0));
  Poss *cont3Poss = new Poss(cont3);
  RealPSet *cont3Set = new RealPSet(cont3Poss);


  Yaxpy *axpy3 = new Yaxpy(SMLAYER, COEFONE);
  axpy3->AddInputs0(2,
		    cont3Set->OutTun(0), q_mnij);
  Poss *axpy3Poss = new Poss(axpy3);
  RealPSet * axpy3Set = new RealPSet(axpy3Poss);

  return axpy3Set;
}

RealPSet* P_jimb_calc(DLANode *r_bmef,
		      DLANode *t_fj, 
		      DLANode *u_jimb, DLANode *w_bmie,
		      DLANode *T_efij, DLANode *x_bmej,
		      DLANode *Tau_efim,
		      const Size big, const Size small)
{
  InputNode *P_jimb = CreateInput4("P_jimb", small, small, small, big);

  Copy *copy1 = new Copy(SMLAYER);
  copy1->AddInputs0(2,
		   u_jimb,
		   P_jimb);
  Poss *copy1Poss = new Poss(copy1);
  RealPSet *copy1Set = new RealPSet(copy1Poss);


  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"bmej","ei","ijmb",(string)"e");
  cont1->AddInputs0(3,
		    x_bmej, t_fj, copy1Set->OutTun(0));
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  
  Contraction *cont2 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"bmie","ej","ijmb",(string)"e");
  cont2->AddInputs0(3,
		    w_bmie, t_fj, cont1Set->OutTun(0));
  Poss *cont2Poss = new Poss(cont2);
  RealPSet *cont2Set = new RealPSet(cont2Poss);

  Contraction *cont4 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"bmef","efij","jimb",(string)"ef");
  cont4->AddInputs0(3,
		    r_bmef, Tau_efim, cont2Set->OutTun(0));
  Poss *cont4Poss = new Poss(cont4);
  RealPSet *cont4Set = new RealPSet(cont4Poss);


  return cont4Set;
}

RealPSet* H_me_calc(DLANode *t_fn, DLANode *v_efmn,
		      const Size big, const Size small)
{
  InputNode *H_me = CreateInput2("H_me", small, big);


  TempVarNode *temp1 = new TempVarNode(v_efmn->DataType(0).GetDist(), "temp1");
  temp1->AddInput(v_efmn,0);
  //  InputNode *temp1 = CreateInput4("temp1", big, big, small, small);

  YAxpPx *axpy1 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "efmn", "efnm");
  axpy1->AddInputs0(3,
		    v_efmn, v_efmn, temp1);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);

  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"efmn","fn","me",(string)"fn");
  cont1->AddInputs0(3,
		    axpy1Set->OutTun(0), t_fn, H_me);
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  return cont1Set;
}

RealPSet* F_ae_calc(DLANode *H_me, DLANode *r_amef,
		    DLANode *t_am,
		    DLANode *v_efmn,
		    DLANode *T_afmn,
		    const Size big, const Size small)
{
  InputNode *F_ae = CreateInput2("F_ae", big, big);


  TempVarNode *temp1 = new TempVarNode(v_efmn->DataType(0).GetDist(), "temp1");
  temp1->AddInput(v_efmn,0);
  TempVarNode *temp2 = new TempVarNode(r_amef->DataType(0).GetDist(), "temp2");
  temp2->AddInput(r_amef,0);

  YAxpPx *axpy1 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "efmn", "efnm");
  axpy1->AddInputs0(3,
		    v_efmn, v_efmn, temp1);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);

  Contraction *cont1 = new Contraction(ABSLAYER,COEFNEGONE,COEFZERO,REAL,"efmn","afmn","ae",(string)"fmn");
  cont1->AddInputs0(3,
		    axpy1Set->OutTun(0), T_afmn, F_ae);
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  YAxpPx *axpy2 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "amef", "amfe");
  axpy2->AddInputs0(3,
		     r_amef, r_amef, temp2);
  Poss *axpy2Poss = new Poss(axpy2);
  RealPSet * axpy2Set = new RealPSet(axpy2Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"amef","fm","ae",(string)"fm");
  cont2->AddInputs0(3,
		    axpy2Set->OutTun(0), t_am, cont1Set->OutTun(0));
  Poss *cont2Poss = new Poss(cont2);
  RealPSet *cont2Set = new RealPSet(cont2Poss);


  Contraction *cont3 = new Contraction(ABSLAYER,COEFNEGONE,COEFONE,REAL,"me","am","ae",(string)"m");
  cont3->AddInputs0(3,
		    H_me, t_am, cont2Set->OutTun(0));
  Poss *cont3Poss = new Poss(cont3);
  RealPSet *cont3Set = new RealPSet(cont3Poss);

  return cont3Set;
}


RealPSet* G_mi_calc(DLANode *H_me, DLANode *u_mnie,
		    DLANode *t_ei,
		    DLANode *v_efmn,
		    DLANode *T_efin,
		    const Size big, const Size small)
{
  InputNode *G_mi = CreateInput2("G_mi", small, small);


  TempVarNode *temp1 = new TempVarNode(v_efmn->DataType(0).GetDist(), "temp1");
  temp1->AddInput(v_efmn,0);
  TempVarNode *temp2 = new TempVarNode(u_mnie->DataType(0).GetDist(), "temp2");
  temp2->AddInput(u_mnie,0);

  YAxpPx *axpy1 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "efmn", "efnm");
  axpy1->AddInputs0(3,
		    v_efmn, v_efmn, temp1);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);

  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"efmn","efin","mi",(string)"efn");
  cont1->AddInputs0(3,
		    axpy1Set->OutTun(0), T_efin, G_mi);
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  YAxpPx *axpy2 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "mnie", "nmie");
  axpy2->AddInputs0(3,
		    u_mnie, u_mnie, temp2);
  Poss *axpy2Poss = new Poss(axpy2);
  RealPSet * axpy2Set = new RealPSet(axpy2Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"mnie","en","mi",(string)"en");
  cont2->AddInputs0(3,
		    axpy2Set->OutTun(0), t_ei, cont1Set->OutTun(0));
  Poss *cont2Poss = new Poss(cont2);
  RealPSet *cont2Set = new RealPSet(cont2Poss);


  Contraction *cont3 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"me","ei","mi",(string)"e");
  cont3->AddInputs0(3,
		    H_me, t_ei, cont2Set->OutTun(0));
  Poss *cont3Poss = new Poss(cont3);
  RealPSet *cont3Set = new RealPSet(cont3Poss);

  return cont3Set;
}

RealPSet* z_ai_calc(DLANode *G_mi,
		    DLANode *H_me, DLANode *U_mnie,
		    DLANode *w_amie,
		    DLANode *x_amei,
		    DLANode *t_am, 
		    DLANode *r_amef, 
		    DLANode *T_aeim,
		    DLANode *Tau_efim,
		    const Size big, const Size small)
{
  InputNode *z_ai = CreateInput2("z_ai", big, small);


  TempVarNode *temp2 = new TempVarNode(r_amef->DataType(0).GetDist(), "temp2");
  temp2->AddInput(r_amef,0);
  TempVarNode *temp3 = new TempVarNode(T_aeim->DataType(0).GetDist(), "temp3");
  temp3->AddInput(T_aeim,0);
  TempVarNode *temp4 = new TempVarNode(w_amie->DataType(0).GetDist(), "temp4");
  temp4->AddInput(w_amie,0);
  TempVarNode *temp5 = new TempVarNode(U_mnie->DataType(0).GetDist(), "temp5");
  temp5->AddInput(U_mnie,0);

  YAxpPx *axpy1 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "amef", "amfe");
  axpy1->AddInputs0(3,
		    r_amef, r_amef, temp2);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"amef","efim","ai",(string)"efm");
  cont2->AddInputs0(3,
		    axpy1Set->OutTun(0), Tau_efim, z_ai);
  Poss *cont2Poss = new Poss(cont2);
  RealPSet *cont2Set = new RealPSet(cont2Poss);


  YAxpPx *axpy2 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "aeim", "aemi");
  axpy2->AddInputs0(3,
		    T_aeim, T_aeim, temp3);
  Poss *axpy2Poss = new Poss(axpy2);
  RealPSet * axpy2Set = new RealPSet(axpy2Poss);

  Contraction *cont3 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"aeim","me","ai",(string)"em");
  cont3->AddInputs0(3,
		    axpy2Set->OutTun(0), H_me, cont2Set->OutTun(0));
  Poss *cont3Poss = new Poss(cont3);
  RealPSet *cont3Set = new RealPSet(cont3Poss);

  YAxpPx *axpy3 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "amei", "amie");
  axpy3->AddInputs0(3,
		    w_amie, x_amei, temp4);
  Poss *axpy3Poss = new Poss(axpy3);
  RealPSet * axpy3Set = new RealPSet(axpy3Poss);

  Contraction *cont4 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"amie","em","ai",(string)"em");
  cont4->AddInputs0(3,
		    axpy3Set->OutTun(0), t_am, cont3Set->OutTun(0));
  Poss *cont4Poss = new Poss(cont4);
  RealPSet *cont4Set = new RealPSet(cont4Poss);


  YAxpPx *axpy4 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "mnie", "nmie");
  axpy4->AddInputs0(3,
		    U_mnie, U_mnie, temp5);
  Poss *axpy4Poss = new Poss(axpy4);
  RealPSet * axpy4Set = new RealPSet(axpy4Poss);

  Contraction *cont5 = new Contraction(ABSLAYER,COEFNEGONE,COEFONE,REAL,"mnie","aemn","ai",(string)"emn");
  cont5->AddInputs0(3,
		    axpy4Set->OutTun(0), T_aeim, cont4Set->OutTun(0));
  Poss *cont5Poss = new Poss(cont5);
  RealPSet *cont5Set = new RealPSet(cont5Poss);

  Contraction *cont6 = new Contraction(ABSLAYER,COEFNEGONE,COEFONE,REAL,"mi","am","ai",(string)"m");
  cont6->AddInputs0(3,
		    G_mi, t_am, cont5Set->OutTun(0));
  Poss *cont6Poss = new Poss(cont6);
  RealPSet *cont6Set = new RealPSet(cont6Poss);

  return cont6Set;
}

RealPSet* Z_abij_calc(DLANode *v_abij, 
		      DLANode *y_abef,
		      DLANode *r_ejab,
		      DLANode *t_am,
		      DLANode *Q_mnij,
		      DLANode *P_ijmb,
		      DLANode *F_ae,
		      DLANode *G_mi, 
		      DLANode *W_bmje,
		      DLANode *X_bmej,
		      DLANode *T_aeim,
		      DLANode *Tau_efim,
		      const Size big, const Size small)
{
  InputNode *Z_abij = CreateInput4("Z_abij", big, big, small, small);


  TempVarNode *temp1 = new TempVarNode(Z_abij->DataType(0).GetDist(), "temp1");
  temp1->AddInput(Z_abij,0);
  TempVarNode *temp2 = new TempVarNode(T_aeim->DataType(0).GetDist(), "temp2");
  temp2->AddInput(T_aeim,0);
  TempVarNode *temp4 = new TempVarNode(T_aeim->DataType(0).GetDist(), "temp4");
  temp4->AddInput(T_aeim,0);
  TempVarNode *accum = new TempVarNode(Z_abij->DataType(0).GetDist(), "accum");
  accum->AddInput(Z_abij,0);



  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"bmej","aemi","abij",(string)"em");
  cont1->AddInputs0(3,
		    X_bmej, T_aeim, temp1);
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  YAxpPx *axpy1 = new YAxpPx(ABSLAYER, COEFNEGONEHALF, COEFNEGONE, "abij", "abji");
  axpy1->AddInputs0(3,
		    cont1Set->OutTun(0), cont1Set->OutTun(0), accum);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);


  YAxpPx *axpy2 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONE, "aeim", "aemi");
  axpy2->AddInputs0(3,
		    T_aeim, T_aeim, temp2);
  Poss *axpy2Poss = new Poss(axpy2);
  RealPSet * axpy2Set = new RealPSet(axpy2Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFONEHALF,COEFONE,REAL,"bmje","aeim","abij",(string)"em");
  cont2->AddInputs0(3,
		    W_bmje, axpy2Set->OutTun(0), axpy1Set->OutTun(0));
  Poss *cont2Poss = new Poss(cont2);
  RealPSet *cont2Set = new RealPSet(cont2Poss);

  Contraction *cont3 = new Contraction(ABSLAYER,COEFNEGONE,COEFONE,REAL,"mi","abmj","abij",(string)"m");
  cont3->AddInputs0(3,
		    G_mi, T_aeim, cont2Set->OutTun(0));
  Poss *cont3Poss = new Poss(cont3);
  RealPSet *cont3Set = new RealPSet(cont3Poss);

  Contraction *cont4 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"ae","ebij","abij",(string)"e");
  cont4->AddInputs0(3,
		    F_ae, T_aeim, cont3Set->OutTun(0));
  Poss *cont4Poss = new Poss(cont4);
  RealPSet *cont4Set = new RealPSet(cont4Poss);

  Contraction *cont5 = new Contraction(ABSLAYER,COEFNEGONE,COEFONE,REAL,"ijmb","am","abij",(string)"m");
  cont5->AddInputs0(3,
		    P_ijmb, t_am, cont4Set->OutTun(0));
  Poss *cont5Poss = new Poss(cont5);
  RealPSet *cont5Set = new RealPSet(cont5Poss);

  Contraction *cont6 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"ejab","ei","abij",(string)"e");
  cont6->AddInputs0(3,
		    r_ejab, t_am, cont5Set->OutTun(0));
  Poss *cont6Poss = new Poss(cont6);
  RealPSet *cont6Set = new RealPSet(cont6Poss);

  YAxpPx *axpy3 = new YAxpPx(ABSLAYER, COEFONE, COEFONE, "abij", "baji");
  axpy3->AddInputs0(3,
		    cont6Set->OutTun(0), cont6Set->OutTun(0), Z_abij);
  Poss *axpy3Poss = new Poss(axpy3);
  RealPSet * axpy3Set = new RealPSet(axpy3Poss);

  Contraction *cont8 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"abef","efij","abij",(string)"ef");
  cont8->AddInputs0(3,
		    y_abef, Tau_efim, axpy3Set->OutTun(0));
  Poss *cont8Poss = new Poss(cont8);
  RealPSet *cont8Set = new RealPSet(cont8Poss);

  Contraction *cont9 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"mnij","abmn","abij",(string)"mn");
  cont9->AddInputs0(3,
		    Q_mnij, Tau_efim, cont8Set->OutTun(0));
  Poss *cont9Poss = new Poss(cont9);
  RealPSet *cont9Set = new RealPSet(cont9Poss);

  Yxpy *axpy4 = new Yxpy(SMLAYER);
  axpy4->AddInputs0(2,
		    v_abij, cont9Set->OutTun(0));
  Poss *axpy4Poss = new Poss(axpy4);
  RealPSet *axpy4Set = new RealPSet(axpy4Poss);

  return axpy4Set;
}

RealPSet* Tau_efmn_calc(DLANode *t_am,
			DLANode *T_aeim,
			const Size big, const Size small)
{
  InputNode *Tau_efmn = CreateInput4("Tau_efmn", big, big, small, small);

  Copy *copy1 = new Copy(SMLAYER);
  copy1->AddInputs0(2,
		    T_aeim,
		    Tau_efmn);
  Poss *copy1Poss = new Poss(copy1);
  RealPSet *copy1Set = new RealPSet(copy1Poss);

  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFONE,REAL,"em","fn","efmn",(string)"");
  cont1->AddInputs0(3,
		    t_am, t_am, copy1Set->OutTun(0));
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  return cont1Set;
}

#endif //DOTENSORS
