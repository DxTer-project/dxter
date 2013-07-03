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
#include "base.h"
#include "comm.h"
#include "distributions.h"

string CommToStr(Comm comm)
{
  switch(comm)
    {
#if DODM
    case(MRCOMM):
      return "g.MRComm()";
    case (MCCOMM):
      return "g.MCComm()";
#elif DOSM||DOSQM
    case(GLOBALCOMM):
      return "GlobalComm";
    case(PROCCOMM):
      return "ProcComm";
    case(L2COMM):
      return "L2Comm";
#endif
    case(CORECOMM):
      return "Seq";
    default:
      throw;
    }
}


// Partial ordering of comm's
// Used to make sure nested parallelization
//  doesn't include an insane setup with
//  L2 parallelization outside of cross-processor
//  parallelization
bool CommAllowedWithin(Comm comm1, Comm comm2)
{
#if DOSM
  switch(comm1)
    {
#if NUMPROCS>1
    case(GLOBALCOMM):
      return comm2 == PROCCOMM || comm2==L2COMM || comm2==CORECOMM;
#endif //NUMPROCS>1

    case(PROCCOMM):
      return comm2==L2COMM || comm2==CORECOMM;

    case(L2COMM):
      return comm2==CORECOMM;

    case(CORECOMM):
      return true;

    default:
      throw;
    }
#else
      throw;
#endif
}
