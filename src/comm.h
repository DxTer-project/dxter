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
#if DOELEM
  MRCOMM,
  MCCOMM,
#elif DOSM||DOSQM
#if NUMPROCS>1
  ALLPROCCOMM,  //Communicate between all cores in system, PROCCOM is subgroup
#if NUML2PERPROC>1
  ALLL2COMM,    //Communicate between all cores in system, L2COMMSUBALLL2 is subgroup
#endif //NUML2PERPROC>1
#endif //NUMPROCS>1
  PROCCOMM,    //Communicate between all cores on processor (e.g. all cores sharing the L3), L2COMM is subgroup
  L2COMM,      //Communicate between all cores sharing my L2, each core is a subgroup
  L2COMMSUBALLL2,    //Communicate between all cores sharing my L2, each core is a subgroup
  L1COMM,
#endif
  CORECOMM,
  BADCOMM
};

Comm MaxComm(Comm comm1, Comm comm2);

string CommToStr(Comm comm);
bool CommAllowedWithin(Comm comm1, Comm comm2);

inline Comm GetSubComm(Comm comm)
{
  switch(comm)
    {
#if DOELEM
      throw;
#elif DOSM||DOSQM
#if NUMPROCS>1
    case(ALLPROCCOMM):
      return PROCCOMM;
#endif //NUMPROCS>1
#if NUML2PERPROC>1
    case (ALLL2COMM):
      return L2COMMSUBALLL2;
#endif //NUML2PERPROC > 1
    case(PROCCOMM):
      return L2COMM;
    case(L2COMM):
      case (L2COMMSUBALLL2):
      return L1COMM;
#endif
    case(CORECOMM):
      throw;
    default:
      throw;
    }
}

inline bool IsImmediateSubComm(Comm commTop, Comm commBott) 
{
  return commBott == GetSubComm(commTop);
}

inline unsigned int NumCoresInComm(Comm comm)
{
  switch(comm)
    {
#if DOELEM
    case(MRCOMM):
      return CVAL;
    case (MCCOMM):
      return RVAL;
#elif DOSM||DOSQM
#if NUMPROCS>1
    case(ALLPROCCOMM):
      return NUMCORESPERL2*NUML2PERPROC*NUMPROCS;
#if NUML2PERPROC>1
    case(ALLL2COMM):
      return NUMCORESPERL2*NUML2PERPROC*NUMPROCS;
#endif //NUML2PERPROC>1
#endif //NUMPROCS>1
    case(PROCCOMM):
      return NUML2PERPROC*NUMCORESPERL2;
    case(L2COMM):
      case(L2COMMSUBALLL2):
      return NUMCORESPERL2;
#endif
    case(CORECOMM):
      return 1;
    default:
      throw;
    }
}

inline unsigned int NumGroupsInComm(Comm comm)
{
  switch(comm)
    {
#if DOSM||DOSQM
#if NUMPROCS>1
    case(ALLPROCCOMM):
      return NUMPROCS;
#if NUML2PERPROC>1
    case(ALLL2COMM):
      return NUMPROCS*NUML2PERPROC;
#endif //NUML2PERPROC>1
#endif //NUMPROCS>1
    case(PROCCOMM):
      return NUML2PERPROC;
    case(L2COMM):
      return NUMCORESPERL2;
#endif
      case(CORECOMM):
        return 1;
    default:
      throw;
    }
}


