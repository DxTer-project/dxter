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
#include "yaxppx.h"
#include "zaxpby.h"
#include "yaxpy.h"
#include "yxpy.h"
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
		      const Size big, const Size small)
{
  InputNode *W_bmje = CreateInput4("W_bmje", big, small, small, big);


  InputNode *temp1 = CreateInput4("temp1", big, big, small, small);
  InputNode *temp2 = CreateInput4("temp2", big, big, small, small);
  InputNode *temp3 = CreateInput4("temp3", small, small, small, big);
  InputNode *temp4 = CreateInput4("temp4", big, small, big, big);

  YAxpPx *axpy0 = new YAxpPx(ABSLAYER, COEFTWO, COEFNEGONEHALF, "bmej", "bmje");
  axpy0->AddInputs0(3,
		     w_bmje, x_bmej, W_bmje);
  Poss *axpy0Poss = new Poss(axpy0);
  RealPSet * axpy0Set = new RealPSet(axpy0Poss);
  

  YAxpPx *axpy1 = new YAxpPx(ABSLAYER, COEFNEGONEHALF, COEFONE, "bfnj", "bfjn");
  axpy1->AddInputs0(3,
		     T_bfnj, T_bfnj, temp1);
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);

  Contraction *cont1 = new Contraction(ABSLAYER,COEFNEGONE,COEFONE,REAL,"bn","fj","bfnj",(string)"");
  cont1->AddInputs0(3,
		    t_fj, t_fj, axpy1Set->OutTun(0));
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  YAxpPx *axpy2 = new YAxpPx(ABSLAYER, COEFNEGONE, COEFTWO, "femn", "fenm");
  axpy2->AddInputs0(3,
		     v_femn, v_femn, temp2);
  Poss *axpy2Poss = new Poss(axpy2);
  RealPSet * axpy2Set = new RealPSet(axpy2Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"femn","bfnj","bmje",(string)"fn");
  cont2->AddInputs0(3,
		    axpy2Set->OutTun(0), cont1Set->OutTun(0), axpy0Set->OutTun(0));
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
		      const Size big, const Size small)
{
  InputNode *X_bmej = CreateInput4("X_bmej", big, small, big, small);


  InputNode *temp1 = CreateInput4("temp1", big, big, small, small);

  Contraction *cont1 = new Contraction(ABSLAYER,COEFONE,COEFZERO,REAL,"bn","fj","bfnj",(string)"");
  cont1->AddInputs0(3,
		    t_fj, t_fj, temp1);
  Poss *cont1Poss = new Poss(cont1);
  RealPSet *cont1Set = new RealPSet(cont1Poss);

  Yaxpy *axpy1 = new Yaxpy(SMLAYER, COEFONEHALF);
  axpy1->AddInputs0(2,
		    T_bfnj, cont1Set->OutTun(0));
  Poss *axpy1Poss = new Poss(axpy1);
  RealPSet * axpy1Set = new RealPSet(axpy1Poss);

  Contraction *cont2 = new Contraction(ABSLAYER,COEFNEGONE,COEFZERO,REAL,"femn","bfnj","bmej",(string)"fn");
  cont2->AddInputs0(3,
		    v_femn, axpy1Set->OutTun(0), X_bmej);
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

  Yxpy *axpy2 = new Yxpy(SMLAYER);
  axpy2->AddInputs0(2,
		    x_bmej, cont4Set->OutTun(0));
  Poss *axpy2Poss = new Poss(axpy2);
  RealPSet * axpy2Set = new RealPSet(axpy2Poss);
  

  return axpy2Set;
}

#endif //DOTENSORS
