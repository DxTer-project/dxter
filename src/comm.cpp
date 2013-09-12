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
#if NUMPROCS>1
    case(ALLPROCCOMM):
      return "GlobalComm";
    case(ALLL2COMM):
      return "AllL2Comm";
#endif //NUML2PERPROC > 1
    case(PROCCOMM):
      return "ProcComm";
    case(L2COMM):
      return "L2Comm";
  case(L1COMM):
    return "L1Comm";
    case(L2COMMSUBALLL2):
      return "L2SubAllL2Comm";
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
    case(ALLPROCCOMM):
      return comm2 == PROCCOMM || comm2==L2COMM || comm2==CORECOMM;
#if NUML2PERPROC>1
    case(ALLL2COMM):
      return comm2==L2COMMSUBALLL2 || comm2==CORECOMM;
#endif //NUML2PERPROC>1
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


Comm MaxComm(Comm comm1, Comm comm2)
{
  switch(comm1)
  {
#if DODM
      throw;
#elif DOSM||DOSQM
#if NUMPROCS>1
    case(ALLPROCCOMM):
      if (comm1 == ALLL2COMM)
        throw;
      else
        return ALLPROCCOMM;
#endif //NUMPROCS>1
#if NUML2PERPROC>1
    case (ALLL2COMM):
      if (comm1 == ALLPROCCOMM)
        throw;
      else
        return ALLL2COMM;
#endif //NUML2PERPROC > 1
    case(PROCCOMM):
      if (comm2 == ALLL2COMM
          || comm2 == ALLPROCCOMM)
        return comm2;
      else
        return PROCCOMM;
    case(L2COMM):
      if (comm2 == ALLL2COMM)
        throw;
      else if (comm2 == ALLPROCCOMM
               || comm2 == PROCCOMM)
        return comm2;
      else
        return L2COMM;
    case (L2COMMSUBALLL2):
      if (comm2 == ALLL2COMM)
        return comm2;
      else if (comm2 == CORECOMM)
        return comm1;
      else if (comm2 == L2COMMSUBALLL2)
        return comm1;
      else if (comm2 == PROCCOMM)
        return CORECOMM; //Let's say we have a Poss with two loops, and internally
                        // those loops have some posses, where one loop has parallelization within
                        // that is PROCCOMM and the other is L2COMMSUBALLL2.
      //When printing the ucrrent Poss, ParallelismWithinCurrent... is called and the max
      // parallelism is determined.  This case can show up, but will never
      // be the best code.  We don't want to throw an except (since we might be printing
      // all codes), so this workaround is horribly nasty so we can print SOMETHING
      // BAM TODO - fix this somehow
      else
        throw;
#endif
    case(CORECOMM):
      return comm2;
    default:
      throw;
  }
}
