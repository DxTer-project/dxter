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
#elsif DOSM
    case(GLOBALCOMM):
      return "GlobalComm";
    case(PROCCOM):
      return "ProcComm";
    case(L2COMM):
      return "L2Comm";
#endif
    default:
      throw;
    }
}

unsigned int NumCoresInComm(Comm comm)
{
  switch(comm)
    {
#if DODM
    case(MRCOMM):
      return CVAL;
    case (MCCOMM):
      return RVAL;
#elsif DOSM
    case(GLOBALCOMM):
      return NUMCORESPERL2*NUML2PERL3*NUML3;
    case(PROCCOM):
      return NUML2PERL3*NUMCORESPERL2;
    case(L2COMM):
      return NUMCORESPERL2;
#endif
    default:
      throw;
    }
}
