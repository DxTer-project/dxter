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
// NOTE: It is possible to simply include "tensormental.hpp" instead
#include "tensormental.hpp"
using namespace tmen;
using namespace std;

#define GRIDORDER 4

template <typename T>
void PrintLocalSizes(const DistTensor<T>& A) 
{
  const Int commRank = mpi::CommRank( mpi::COMM_WORLD );
  if (commRank == 0) {
    for (Unsigned i = 0; i < A.Order(); ++i) {
      cout << i << " is " << A.LocalDimension(i) << endl;
    }
  }
}


template <typename T>
void GatherAllModes(const DistTensor<T>& A, DistTensor<T>& B)
{
  DistTensor<T> *tmp = NULL;
  //  DistTensor<T> *tmp = new DistTensor<T>(A.TensorDist(), A.Grid());
  //  *tmp = A;
  
  const TensorDistribution dist = A.TensorDist();

  for (Unsigned mode = 0; mode < A.Order(); ++mode) {
    ModeDistribution modeDist = dist[mode];
    if (!(modeDist.empty())) {
      TensorDistribution newDist = (tmp ? tmp->TensorDist() : A.TensorDist());
      ModeDistribution newIgnoreDist = newDist[newDist.size() - 1];
      newIgnoreDist.insert(newIgnoreDist.end(), modeDist.begin(), modeDist.end());
      newDist[newDist.size() - 1] = newIgnoreDist;
      modeDist.clear();
      newDist[mode] = modeDist;
      DistTensor<T> *tmp2 = new DistTensor<T>(newDist, A.Grid());
      if (!tmp) {
	tmp2->GatherToOneRedistFrom(A, mode);
      }
      else {
	tmp2->GatherToOneRedistFrom(*tmp, mode);
	delete tmp;
      }
      tmp = tmp2;
    }
  }

  if (tmp) {
    
    if (TensorDistToString(B.TensorDist()) != TensorDistToString(tmp->TensorDist())) {
      cout << TensorDistToString(B.TensorDist()) << endl;
      cout << TensorDistToString(tmp->TensorDist()) << endl;
      throw;
    }

    B = *tmp;
    delete tmp;
  }
  else {
    B = A;
  }
}


void Usage(){
  std::cout << "./DistTensor <gridDim0> <gridDim1> ... \n";
    std::cout << "<gridDimK>   : dimension of mode-K of grid\n";
}

typedef struct Arguments{
  ObjShape gridShape;
  Unsigned nProcs;
} Params;

void ProcessInput(int argc,  char** const argv, Params& args){
    Unsigned i;
    Unsigned argCount = 0;
    if(argCount + 1 >= argc){
        std::cerr << "Missing required gridOrder argument\n";
        Usage();
        throw ArgException();
    }

    if(argCount + GRIDORDER >= argc){
        std::cerr << "Missing required grid dimensions\n";
        Usage();
        throw ArgException();
    }

    args.gridShape.resize(GRIDORDER);
    args.nProcs = 1;
    for(int i = 0; i < GRIDORDER; i++){
        int gridDim = atoi(argv[++argCount]);
        if(gridDim <= 0){
            std::cerr << "Grid dim must be greater than 0\n";
            Usage();
            throw ArgException();
        }
	args.nProcs *= gridDim;
        args.gridShape[i] = gridDim;
    }
}

template<typename T>
void
Set(DistTensor<T>& A)
{
  Unsigned order = A.Order();
  Location loc(order);
  std::fill(loc.begin(), loc.end(), 0);
  Unsigned ptr = 0;
  bool stop = false;

  while(!stop){
    A.Set(loc, rand());
    if (loc.size() == 0)
      break;

    //Update
    loc[ptr]++;
    while(loc[ptr] == A.Dimension(ptr)){
      loc[ptr] = 0;
      ptr++;
      if(ptr == order){
	stop = true;
	break;
      }else{
	loc[ptr]++;
      }
    }
    ptr = 0;
  }
}

template<typename T>
void
DistTensorTest( const Grid& g )
{
#ifndef RELEASE
  CallStackEntry entry("DistTensorTest");
#endif
  Unsigned i;
  const Int commRank = mpi::CommRank( mpi::COMM_WORLD );
  const Unsigned gridOrder = 4;

  ObjShape tempShape;
  TensorDistribution dist____N_D_0_1_2_3 = tmen::StringToTensorDist("[]|(0,1,2,3)");
  TensorDistribution dist__D_0__D_1__D_2__D_3 = tmen::StringToTensorDist("[(0),(1),(2),(3)]");
  TensorDistribution dist__D_0__D_1__D_3__D_2 = tmen::StringToTensorDist("[(0),(1),(3),(2)]");
  //E_MP2[D0,D1,D2,D3]
  DistTensor<double> E_MP2__D_0__D_1__D_2__D_3( dist__D_0__D_1__D_2__D_3, g );
  //E_MP2[] | {0,1,2,3}
  DistTensor<double> E_MP2____N_D_0_1_2_3( dist____N_D_0_1_2_3, g );
  //axppx1_temp[D0,D1,D2,D3]
  DistTensor<double> axppx1_temp__D_0__D_1__D_2__D_3( dist__D_0__D_1__D_2__D_3, g );
  //t_efmn[D0,D1,D2,D3]
  DistTensor<double> t_efmn__D_0__D_1__D_2__D_3( dist__D_0__D_1__D_2__D_3, g );
  //v_efmn[D0,D1,D2,D3]
  DistTensor<double> v_efmn__D_0__D_1__D_2__D_3( dist__D_0__D_1__D_2__D_3, g );
  //v_efmn[D0,D1,D3,D2]
  DistTensor<double> v_efmn__D_0__D_1__D_3__D_2( dist__D_0__D_1__D_3__D_2, g );
  ModeArray modes_0_1_2_3;
  modes_0_1_2_3.push_back(0);
  modes_0_1_2_3.push_back(1);
  modes_0_1_2_3.push_back(2);
  modes_0_1_2_3.push_back(3);
  ModeArray modes_2;
  modes_2.push_back(2);
  ModeArray modes_2_3;
  modes_2_3.push_back(2);
  modes_2_3.push_back(3);
  ModeArray modes_3;
  modes_3.push_back(3);
  ModeArray modes_3_2;
  modes_3_2.push_back(3);
  modes_3_2.push_back(2);
  IndexArray indices_efmn( 4 );
  indices_efmn[0] = 'e';
  indices_efmn[1] = 'f';
  indices_efmn[2] = 'm';
  indices_efmn[3] = 'n';
  std::vector<ModeArray> modeArrayArray___2___3;
  modeArrayArray___2___3.push_back(modes_2);
  modeArrayArray___2___3.push_back(modes_3);
  Permutation perm_0_1_3_2;
  perm_0_1_3_2.push_back(0);
  perm_0_1_3_2.push_back(1);
  perm_0_1_3_2.push_back(3);
  perm_0_1_3_2.push_back(2);
  // v_efmn has 4 dims
  //	Starting distribution: [D0,D1,D2,D3] or _D_0__D_1__D_2__D_3
  ObjShape v_efmn__D_0__D_1__D_2__D_3_tempShape;
  v_efmn__D_0__D_1__D_2__D_3_tempShape.push_back( 53 );
  v_efmn__D_0__D_1__D_2__D_3_tempShape.push_back( 53 );
  v_efmn__D_0__D_1__D_2__D_3_tempShape.push_back( 5 );
  v_efmn__D_0__D_1__D_2__D_3_tempShape.push_back( 5 );
  v_efmn__D_0__D_1__D_2__D_3.ResizeTo( v_efmn__D_0__D_1__D_2__D_3_tempShape );
  MakeUniform( v_efmn__D_0__D_1__D_2__D_3 );
  DistTensor<T> v_efmn_local( tmen::StringToTensorDist("[(),(),(),()]|(0,1,2,3)"), g );
  GatherAllModes( v_efmn__D_0__D_1__D_2__D_3, v_efmn_local );
  // axppx1_temp has 4 dims
  //	Starting distribution: [D0,D1,D2,D3] or _D_0__D_1__D_2__D_3
  ObjShape axppx1_temp__D_0__D_1__D_2__D_3_tempShape;
  axppx1_temp__D_0__D_1__D_2__D_3_tempShape.push_back( 53 );
  axppx1_temp__D_0__D_1__D_2__D_3_tempShape.push_back( 53 );
  axppx1_temp__D_0__D_1__D_2__D_3_tempShape.push_back( 5 );
  axppx1_temp__D_0__D_1__D_2__D_3_tempShape.push_back( 5 );
  axppx1_temp__D_0__D_1__D_2__D_3.ResizeTo( axppx1_temp__D_0__D_1__D_2__D_3_tempShape );
  MakeUniform( axppx1_temp__D_0__D_1__D_2__D_3 );
  DistTensor<T> axppx1_temp_local( tmen::StringToTensorDist("[(),(),(),()]|(0,1,2,3)"), g );
  GatherAllModes( axppx1_temp__D_0__D_1__D_2__D_3, axppx1_temp_local );
  // t_efmn has 4 dims
  //	Starting distribution: [D0,D1,D2,D3] or _D_0__D_1__D_2__D_3
  ObjShape t_efmn__D_0__D_1__D_2__D_3_tempShape;
  t_efmn__D_0__D_1__D_2__D_3_tempShape.push_back( 53 );
  t_efmn__D_0__D_1__D_2__D_3_tempShape.push_back( 53 );
  t_efmn__D_0__D_1__D_2__D_3_tempShape.push_back( 5 );
  t_efmn__D_0__D_1__D_2__D_3_tempShape.push_back( 5 );
  t_efmn__D_0__D_1__D_2__D_3.ResizeTo( t_efmn__D_0__D_1__D_2__D_3_tempShape );
  MakeUniform( t_efmn__D_0__D_1__D_2__D_3 );
  DistTensor<T> t_efmn_local( tmen::StringToTensorDist("[(),(),(),()]|(0,1,2,3)"), g );
  GatherAllModes( t_efmn__D_0__D_1__D_2__D_3, t_efmn_local );
  // scalar input has 0 dims
  //	Starting distribution: [] | {0,1,2,3} or ___N_D_0_1_2_3
  ObjShape E_MP2____N_D_0_1_2_3_tempShape;
  E_MP2____N_D_0_1_2_3.ResizeTo( E_MP2____N_D_0_1_2_3_tempShape );
  MakeUniform( E_MP2____N_D_0_1_2_3 );
  DistTensor<T> E_MP2_local( tmen::StringToTensorDist("[]|(0,1,2,3)"), g );
  GatherAllModes( E_MP2____N_D_0_1_2_3, E_MP2_local );
  //**** (out of 1)
  //------------------------------------//

  // v_efmn[D0,D1,D3,D2] <- v_efmn[D0,D1,D2,D3]
  v_efmn__D_0__D_1__D_3__D_2.AllToAllRedistFrom( v_efmn__D_0__D_1__D_2__D_3, modes_2_3, modes_3_2, modeArrayArray___2___3 );
  YAxpPx( 2.0, v_efmn__D_0__D_1__D_2__D_3, -1.0, v_efmn__D_0__D_1__D_3__D_2, perm_0_1_3_2, axppx1_temp__D_0__D_1__D_2__D_3 );

  //------------------------------------//

  //****
  //**** (out of 1)
  //------------------------------------//

  tempShape = E_MP2____N_D_0_1_2_3.Shape();
  tempShape.push_back( g.Shape()[0] );
  tempShape.push_back( g.Shape()[1] );
  tempShape.push_back( g.Shape()[2] );
  tempShape.push_back( g.Shape()[3] );
  E_MP2__D_0__D_1__D_2__D_3.ResizeTo( tempShape );
  // 1.0 * axppx1_temp[D0,D1,D2,D3]_efmn * t_efmn[D0,D1,D2,D3]_efmn + 0.0 * E_MP2[D0,D1,D2,D3]_efmn
  LocalContract(1.0, axppx1_temp__D_0__D_1__D_2__D_3.LockedTensor(), indices_efmn,
		t_efmn__D_0__D_1__D_2__D_3.LockedTensor(), indices_efmn,
		0.0, E_MP2__D_0__D_1__D_2__D_3.Tensor(), indices_efmn);
  // E_MP2[] | {0,1,2,3} <- E_MP2[D0,D1,D2,D3] (with SumScatter on (D0)(D1)(D2)(D3))
  E_MP2____N_D_0_1_2_3.ReduceToOneRedistFrom( E_MP2__D_0__D_1__D_2__D_3, modes_0_1_2_3 );

  //------------------------------------//

//****
/*
      DistTensor<T> diffTensor( tmen::StringToTensorDist("[(),(),()]|(0,1,2,3)"), g );    
      diffTensor.ResizeTo(C_local);
      Diff( C_local.LockedTensor(), C_local_comparison.LockedTensor(), diffTensor.Tensor() );

      if (commRank == 0) {
      cout << "Norm of distributed is " << Norm(C_local_comparison.LockedTensor()) << endl;
      cout << "Norm of local is " << Norm(C_local.LockedTensor()) << endl;
      cout << "Norm is " << Norm(diffTensor.LockedTensor()) << endl;
      }

    */
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    Unsigned i;
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::CommRank( comm );
    const Int commSize = mpi::CommSize( comm );
    printf("My Rank: %d\n", commRank);
    try
    {
        Params args;

        ProcessInput(argc, argv, args);

        if(commRank == 0 && commSize != args.nProcs){
            std::cerr << "program not started with correct number of processes\n";
	    std::cerr << commSize << " vs " << args.nProcs << std::endl;
            Usage();
            throw ArgException();
        }

        if(commRank == 0){
            printf("Creating %d", args.gridShape[0]);
            for(i = 1; i < GRIDORDER; i++)
                printf(" x %d", args.gridShape[i]);
            printf(" grid\n");
        }

        const Grid g( comm, args.gridShape );

        if( commRank == 0 )
        {
            std::cout << "------------------" << std::endl
                      << "Testing with doubles:" << std::endl
                      << "------------------" << std::endl;
        }
        DistTensorTest<double>( g );

    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    //printf("Completed\n");
    return 0;
}
