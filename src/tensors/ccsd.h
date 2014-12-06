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
		      const Size big, const Size small);

RealPSet* X_bmej_calc(DLANode *x_bmej, DLANode *r_bmef,
		      DLANode *t_fj, 
		      DLANode *u_mnje, DLANode *v_femn,
		      DLANode *T_bfnj,
		      const Size big, const Size small);

RealPSet* U_mnie_calc(DLANode *t_fj, 
		      DLANode *u_mnie, DLANode *v_femn,
		      const Size big, const Size small);

RealPSet* Q_mnij_calc(DLANode *q_mnij,
		      DLANode *t_fj, 
		      DLANode *u_mnie, DLANode *v_femn,
		      DLANode *T_bfnj,
		      const Size big, const Size small);

RealPSet* P_jimb_calc(DLANode *r_bmef,
		      DLANode *t_fj, 
		      DLANode *u_jimb, DLANode *w_bmie,
		      DLANode *T_efij, DLANode *x_bmej,
		      const Size big, const Size small);

RealPSet* H_me_calc(DLANode *t_fn, DLANode *v_efmn,
		      const Size big, const Size small);



InputNode *CreateInput2(string name, Size size1, Size size2);
InputNode *CreateInput4(string name, Size size1, Size size2, Size size3, Size size4);


#endif //DOTENSORS
