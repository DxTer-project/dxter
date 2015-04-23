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

#include "driverMenu.h"

#if DOLLDLA

void PrintMainMenu()
{
  cout <<"\n\nWelcome to DxTer LLDLA! Please choose an option below:\n\n";
  cout << "./driver arg1 arg2 ...\n";
  cout <<"\n";
  cout <<"arg1 == 0   -> View benchmarks\n";
  cout <<"\n";
  cout <<"Single Operation Examples\n";
  cout <<"         3  -> Dot prod F/D M\n";
  cout <<"         4  -> Matrix add F/D M N\n";
  cout <<"         5  -> Matrix vector multiply N/T F/D M N\n";
  cout <<"         6  -> Scalar vector multiply C/R F/D M\n";
  cout <<"         7  -> Vector matrix multiply F/D M N\n";
  cout <<"         8  -> Scalar matrix multiply F/D M N\n";
  cout <<"         9  -> Vector add C/R F/D M\n";
  cout <<"        15  -> Gen Size Col Vector SVMul F/D M\n";
  cout <<"        33  -> MMMul F/D M N P\n";
  cout <<"\n";
  cout <<"BLAS Examples\n";
  cout <<"         1  -> Gemm  N/T N/T F/D M N P\n";
  cout <<"        14  -> Gemv N/T F/D M N\n";
  cout <<"        16  -> Axpy C/R F/D M\n";
  cout <<"\n";
  cout <<"Miscellaneous Examples\n";
  cout <<"         2  -> Double Gemm  N/T N/T F/D M N P K\n";
  cout <<"        10  -> Vector add twice F/D M\n";
  cout <<"        11  -> Vector matrix vector multiply F/D M N\n";
  cout <<"        12  -> Matrix add twice F/D M N\n";
  cout <<"        13  -> Matrix vector multiply twice F/D M N P\n";
  cout <<"        17  -> alpha*(A0 + A1)^T*B + beta*C^T F/D M N P\n";
  cout <<"        18  -> alpha*A*x + beta*B*x F/D M N\n";
  cout <<"        31  -> y = alpha*x + beta*(z + y) F/D M\n";
  cout <<"        32  -> y = (A + B^T)*x F/D M N\n";
  cout <<"        34  -> q = A*p; s = A^T*r F/D M N\n";
  cout <<"\n";
  cout <<"Node test examples\n";
  cout <<"        19  -> Set to zero test (y <- Ax) F/D M N\n";
  cout <<"        20  -> Pack test F/D M\n";
  cout <<"        21  -> Copy test F/D M N\n";
  cout <<"        22  -> Vertical partition recombine test F/D M\n";
  cout <<"        23  -> Horizontal partition recombine test F/D M\n";
  cout <<"        24  -> Vertical refined pack test F/D M\n";
  cout <<"        25  -> Vertical pack unpack test F/D M\n";
  cout <<"        26  -> 2D vertical pack unpack test F/D M N\n";
  cout <<"        27  -> 2D vertical unpack test F/D M N\n";
  cout <<"        28  -> 2D horizontal unpack test F/D M N\n";
  cout <<"        29  -> 2D horizontal copy test F/D M N\n";
  cout <<"\n";
  cout <<"Automated tests\n";
  cout <<"        30  -> Basic examples, no runtime evaluation\n";
  cout <<"\n";
}

#endif // DOLLDLA
