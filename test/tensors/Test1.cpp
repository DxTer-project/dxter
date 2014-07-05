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
    Unsigned counter = 0;
    bool stop = false;

    while(!stop){
        A.Set(loc, counter);
	if (loc.size() == 0)
	  break;

        //Update
        counter++;
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

    ObjShape shapes3(3, 10);
    ObjShape shapes4(4, 10);





ObjShape tempShape;
TensorDistribution dist__S__D_0__D_3 = tmen::StringToTensorDist("[(),(0),(3)]");
TensorDistribution dist__S__D_1_2__S__D_0__D_3 = tmen::StringToTensorDist("[(),(1,2),(),(0),(3)]");
TensorDistribution dist__S__D_2_0__D_3 = tmen::StringToTensorDist("[(),(2,0),(3)]");
TensorDistribution dist__S__D_2_1__S__D_0__D_3 = tmen::StringToTensorDist("[(),(2,1),(),(0),(3)]");
TensorDistribution dist__S__D_1__S__D_0__D_3 = tmen::StringToTensorDist("[(),(1),(),(0),(3)]");
TensorDistribution dist__S__D_1__D_2__D_0__D_3 = tmen::StringToTensorDist("[(),(1),(2),(0),(3)]");
TensorDistribution dist__S__D_2__S__D_0__D_3 = tmen::StringToTensorDist("[(),(2),(),(0),(3)]");
TensorDistribution dist__S__D_2__D_3 = tmen::StringToTensorDist("[(),(2),(3)]");
TensorDistribution dist__S__D_0_2__D_3 = tmen::StringToTensorDist("[(),(0,2),(3)]");
TensorDistribution dist__D_0__D_1__D_2__D_3 = tmen::StringToTensorDist("[(0),(1),(2),(3)]");
TensorDistribution dist__D_0__D_2__D_3 = tmen::StringToTensorDist("[(0),(2),(3)]");
TensorDistribution dist__D_1__D_2__S__D_0__D_3 = tmen::StringToTensorDist("[(1),(2),(),(0),(3)]");
TensorDistribution dist__D_0_1__D_2__D_3 = tmen::StringToTensorDist("[(0,1),(2),(3)]");
TensorDistribution dist__D_1_0__D_2__S__D_3 = tmen::StringToTensorDist("[(1,0),(2),(),(3)]");
TensorDistribution dist__D_1_0__D_2__D_3 = tmen::StringToTensorDist("[(1,0),(2),(3)]");
	//A[D01,D2,D3]
DistTensor<double> A__D_0_1__D_2__D_3( dist__D_0_1__D_2__D_3, g );
	//A[D0,D2,D3]
DistTensor<double> A__D_0__D_2__D_3( dist__D_0__D_2__D_3, g );
	//A[*,D02,D3]
DistTensor<double> A__S__D_0_2__D_3( dist__S__D_0_2__D_3, g );
	//A[*,D0,D3]
DistTensor<double> A__S__D_0__D_3( dist__S__D_0__D_3, g );
	//A[*,D20,D3]
DistTensor<double> A__S__D_2_0__D_3( dist__S__D_2_0__D_3, g );
	//A[*,D2,D3]
DistTensor<double> A__S__D_2__D_3( dist__S__D_2__D_3, g );
	//B[D0,D1,D2,D3]
DistTensor<double> B__D_0__D_1__D_2__D_3( dist__D_0__D_1__D_2__D_3, g );
	//C[D01,D2,D3]
DistTensor<double> C__D_0_1__D_2__D_3( dist__D_0_1__D_2__D_3, g );
	//C[D10,D2,D3]
DistTensor<double> C__D_1_0__D_2__D_3( dist__D_1_0__D_2__D_3, g );
	//C[D10,D2,*,D3]
DistTensor<double> C__D_1_0__D_2__S__D_3( dist__D_1_0__D_2__S__D_3, g );
	//C[D1,D2,*,D0,D3]
DistTensor<double> C__D_1__D_2__S__D_0__D_3( dist__D_1__D_2__S__D_0__D_3, g );
	//C[*,D12,*,D0,D3]
DistTensor<double> C__S__D_1_2__S__D_0__D_3( dist__S__D_1_2__S__D_0__D_3, g );
	//C[*,D1,D2,D0,D3]
DistTensor<double> C__S__D_1__D_2__D_0__D_3( dist__S__D_1__D_2__D_0__D_3, g );
	//C[*,D1,*,D0,D3]
DistTensor<double> C__S__D_1__S__D_0__D_3( dist__S__D_1__S__D_0__D_3, g );
	//C[*,D21,*,D0,D3]
DistTensor<double> C__S__D_2_1__S__D_0__D_3( dist__S__D_2_1__S__D_0__D_3, g );
	//C[*,D2,*,D0,D3]
DistTensor<double> C__S__D_2__S__D_0__D_3( dist__S__D_2__S__D_0__D_3, g );
ModeArray modes_0;
modes_0.push_back(0);
ModeArray modes_0_1;
modes_0_1.push_back(0);
modes_0_1.push_back(1);
ModeArray modes_1;
modes_1.push_back(1);
ModeArray modes_1_0;
modes_1_0.push_back(1);
modes_1_0.push_back(0);
ModeArray modes_1_2;
modes_1_2.push_back(1);
modes_1_2.push_back(2);
ModeArray modes_2;
modes_2.push_back(2);
ModeArray modes_2_0;
modes_2_0.push_back(2);
modes_2_0.push_back(0);
IndexArray indices_acd( 3 );
indices_acd[0] = 'a';
indices_acd[1] = 'c';
indices_acd[2] = 'd';
IndexArray indices_aefcd( 5 );
indices_aefcd[0] = 'a';
indices_aefcd[1] = 'e';
indices_aefcd[2] = 'f';
indices_aefcd[3] = 'c';
indices_aefcd[4] = 'd';
IndexArray indices_cefd( 4 );
indices_cefd[0] = 'c';
indices_cefd[1] = 'e';
indices_cefd[2] = 'f';
indices_cefd[3] = 'd';
// A input has 3 dims
//	Starting distribution: [D01,D2,D3] or _D_0_1__D_2__D_3
ObjShape A__D_0_1__D_2__D_3_tempShape;
A__D_0_1__D_2__D_3_tempShape.push_back( 10 );
A__D_0_1__D_2__D_3_tempShape.push_back( 10 );
A__D_0_1__D_2__D_3_tempShape.push_back( 10 );
A__D_0_1__D_2__D_3.ResizeTo( A__D_0_1__D_2__D_3_tempShape );
Set( A__D_0_1__D_2__D_3 );
DistTensor<T> A_local( tmen::StringToTensorDist("[(),(),()]|(0,1,2,3)"), g );
GatherAllModes( A__D_0_1__D_2__D_3, A_local );
// B input has 4 dims
//	Starting distribution: [D0,D1,D2,D3] or _D_0__D_1__D_2__D_3
ObjShape B__D_0__D_1__D_2__D_3_tempShape;
B__D_0__D_1__D_2__D_3_tempShape.push_back( 10 );
B__D_0__D_1__D_2__D_3_tempShape.push_back( 10 );
B__D_0__D_1__D_2__D_3_tempShape.push_back( 10 );
B__D_0__D_1__D_2__D_3_tempShape.push_back( 10 );
B__D_0__D_1__D_2__D_3.ResizeTo( B__D_0__D_1__D_2__D_3_tempShape );
Set( B__D_0__D_1__D_2__D_3 );
DistTensor<T> B_local( tmen::StringToTensorDist("[(),(),(),()]|(0,1,2,3)"), g );
GatherAllModes( B__D_0__D_1__D_2__D_3, B_local );
// C input has 3 dims
//	Starting distribution: [D01,D2,D3] or _D_0_1__D_2__D_3
ObjShape C__D_0_1__D_2__D_3_tempShape;
C__D_0_1__D_2__D_3_tempShape.push_back( 10 );
C__D_0_1__D_2__D_3_tempShape.push_back( 10 );
C__D_0_1__D_2__D_3_tempShape.push_back( 10 );
C__D_0_1__D_2__D_3.ResizeTo( C__D_0_1__D_2__D_3_tempShape );
Set( C__D_0_1__D_2__D_3 );
DistTensor<T> C_local( tmen::StringToTensorDist("[(),(),()]|(0,1,2,3)"), g );
GatherAllModes( C__D_0_1__D_2__D_3, C_local );
//**** (out of 4)
	//------------------------------------//

	tempShape = C__D_0_1__D_2__D_3.Shape();
	tempShape.push_back( g.Shape()[0] );
	tempShape.push_back( g.Shape()[3] );
	C__S__D_1__D_2__D_0__D_3.ResizeTo( tempShape );
	   // C[D10,D2,D3] <- C[D01,D2,D3]
	C__D_1_0__D_2__D_3.PermutationRedistFrom( C__D_0_1__D_2__D_3, 0, modes_0_1 );
	tempShape = C__D_1_0__D_2__D_3.Shape();
	tempShape.push_back( g.Shape()[3] );
	C__D_1_0__D_2__S__D_3.ResizeTo( tempShape );
	//**** (out of 2)
		//------------------------------------//

		   // A[D0,D2,D3] <- A[D01,D2,D3]
		A__D_0__D_2__D_3.AllGatherRedistFrom( A__D_0_1__D_2__D_3, 0, modes_1 );
		   // A[*,D2,D3] <- A[D0,D2,D3]
		A__S__D_2__D_3.AllGatherRedistFrom( A__D_0__D_2__D_3, 0, modes_0 );
		   // A[*,D20,D3] <- A[*,D2,D3]
		A__S__D_2_0__D_3.LocalRedistFrom( A__S__D_2__D_3, 1, modes_0 );
		   // A[*,D02,D3] <- A[*,D20,D3]
		A__S__D_0_2__D_3.PermutationRedistFrom( A__S__D_2_0__D_3, 1, modes_2_0 );
		   // A[*,D0,D3] <- A[*,D02,D3]
		A__S__D_0__D_3.AllGatherRedistFrom( A__S__D_0_2__D_3, 1, modes_2 );

		//------------------------------------//

	//****
	   // 1.0 * A[*,D0,D3]_acd * B[D0,D1,D2,D3]_cefd + 0.0 * C[*,D1,D2,D0,D3]_aefcd
	LocalContract(1.0, A__S__D_0__D_3.LockedTensor(), indices_acd,
		B__D_0__D_1__D_2__D_3.LockedTensor(), indices_cefd,
		0.0, C__S__D_1__D_2__D_0__D_3.Tensor(), indices_aefcd);
	   // C[*,D1,*,D0,D3] <- C[*,D1,D2,D0,D3]
	C__S__D_1__S__D_0__D_3.AllGatherRedistFrom( C__S__D_1__D_2__D_0__D_3, 2, modes_2 );
	   // C[*,D12,*,D0,D3] <- C[*,D1,*,D0,D3]
	C__S__D_1_2__S__D_0__D_3.LocalRedistFrom( C__S__D_1__S__D_0__D_3, 1, modes_2 );
	   // C[*,D21,*,D0,D3] <- C[*,D12,*,D0,D3]
	C__S__D_2_1__S__D_0__D_3.PermutationRedistFrom( C__S__D_1_2__S__D_0__D_3, 1, modes_1_2 );
	   // C[*,D2,*,D0,D3] <- C[*,D21,*,D0,D3]
	C__S__D_2__S__D_0__D_3.AllGatherRedistFrom( C__S__D_2_1__S__D_0__D_3, 1, modes_1 );
	   // C[D1,D2,*,D0,D3] <- C[*,D2,*,D0,D3]
	C__D_1__D_2__S__D_0__D_3.LocalRedistFrom( C__S__D_2__S__D_0__D_3, 0, modes_1 );
	   // C[D10,D2,*,D3] <- C[D1,D2,*,D0,D3] (with SumScatter on D0)
	C__D_1_0__D_2__S__D_3.ReduceScatterRedistFrom( C__D_1__D_2__S__D_0__D_3, 3, 0 );
	   // C[D10,D2,D3] <- C[D10,D2,*,D3] (with SumScatter on D3)
	C__D_1_0__D_2__D_3.ReduceScatterUpdateRedistFrom( C__D_1_0__D_2__S__D_3, 1.0, 3, 2 );
	   // C[D01,D2,D3] <- C[D10,D2,D3]
	C__D_0_1__D_2__D_3.PermutationRedistFrom( C__D_1_0__D_2__D_3, 0, modes_1_0 );



    IndexArray indices_aef( 3 );
    indices_aef[0] = 'a';
    indices_aef[1] = 'e';
    indices_aef[2] = 'f';

    LocalContractAndLocalEliminate(1.0, A_local.LockedTensor(), indices_acd,
			      B_local.LockedTensor(), indices_cefd,
			      1.0, C_local.Tensor(), indices_aef);

    DistTensor<T> C_local_comparison( tmen::StringToTensorDist("[(),(),()]|(0,1,2,3)"), g );    
    GatherAllModes(C__D_0_1__D_2__D_3, C_local_comparison);

    DistTensor<T> diffTensor( tmen::StringToTensorDist("[(),(),()]|(0,1,2,3)"), g );    
    diffTensor.ResizeTo(C_local);
    Diff( C_local.LockedTensor(), C_local_comparison.LockedTensor(), diffTensor.Tensor() );

    if (commRank == 0) {
      cout << "Norm of distributed is " << Norm(C_local_comparison.LockedTensor()) << endl;
      cout << "Norm of local is " << Norm(C_local.LockedTensor()) << endl;
      cout << "Norm is " << Norm(diffTensor.LockedTensor()) << endl;
    }
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
