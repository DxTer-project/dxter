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

#pragma once

#include "layers.h"

#if DOTENSORS

#include "transform.h"
#include "DLAOp.h"
#include "tensorRedist.h"

RealPSet* W_bmje_calc(DLANode *w_bmje, DLANode *x_bmej,
		      DLANode *r_bmfe, DLANode *t_fj,
		      DLANode *u_mnje, DLANode *v_femn,
		      DLANode *T_bfnj,
		      DLANode *Tau_efim,
		      const Size big, const Size small);

RealPSet* X_bmej_calc(DLANode *x_bmej, DLANode *r_bmef,
		      DLANode *t_fj, 
		      DLANode *u_mnje, DLANode *v_femn,
		      DLANode *T_bfnj,
		      DLANode *Tau_efim,
		      const Size big, const Size small);

RealPSet* U_mnie_calc(DLANode *t_fj, 
		      DLANode *u_mnie, DLANode *v_femn,
		      const Size big, const Size small);

RealPSet* Q_mnij_calc(DLANode *q_mnij,
		      DLANode *t_fj, 
		      DLANode *u_mnie, DLANode *v_femn,
		      DLANode *T_bfnj,
		      DLANode *Tau_efim,
		      const Size big, const Size small);

RealPSet* P_jimb_calc(DLANode *r_bmef,
		      DLANode *t_fj, 
		      DLANode *u_jimb, DLANode *w_bmie,
		      DLANode *T_efij, DLANode *x_bmej,
		      DLANode *Tau_efim,
		      const Size big, const Size small);

RealPSet* H_me_calc(DLANode *t_fn, DLANode *v_efmn,
		      const Size big, const Size small);

RealPSet* F_ae_calc(DLANode *H_me, DLANode *r_amef,
		    DLANode *t_am,
		    DLANode *v_efmn,
		    DLANode *T_afmn,
		    const Size big, const Size small);

RealPSet* G_mi_calc(DLANode *H_me, DLANode *u_mnie,
		    DLANode *t_ei,
		    DLANode *v_efmn,
		    DLANode *T_efin,
		    const Size big, const Size small);


RealPSet* z_ai_calc(DLANode *G_mi,
		    DLANode *H_me, DLANode *U_mnie,
		    DLANode *w_amie,
		    DLANode *x_amei,
		    DLANode *t_am, 
		    DLANode *r_amef, 
		    DLANode *T_aeim,
		    DLANode *Tau_efim,
		    const Size big, const Size small);


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
		      const Size big, const Size small);

RealPSet* Tau_efmn_calc(DLANode *t_am,
			DLANode *T_aeim,
			const Size big, const Size small);



InputNode *CreateInput2(string name, Size size1, Size size2);
InputNode *CreateInput4(string name, Size size1, Size size2, Size size3, Size size4);
InputNode *CreateInput6(string name, Size size1, Size size2, Size size3, Size size4, Size size5, Size size6);
InputNode *CreateInput3(string name, Size size1, Size size2, Size size3);


#endif //DOTENSORS
