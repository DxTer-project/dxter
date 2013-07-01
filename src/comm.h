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



#pragma once

//#include "transform.h"
//#include "DLANode.h"

enum Comm {
#if DODM
  MRCOMM,
  MCCOMM,
#elif DOSM||DOSQM
  GLOBALCOMM,  //Communicate between all cores in system
  PROCCOMM,    //Communicate between all cores on processor (e.g. all cores sharing the L3)
  L2COMM,      //Communicate between all cores sharing my L2
#endif
  CORECOMM,
  BADCOMM
};

string CommToStr(Comm comm);
unsigned int NumCoresInComm(Comm comm);
unsigned int NumGroupsInComm(Comm comm);

bool CommGroupGreaterThan(Comm comm1, Comm comm2);
