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
void GatherAllModes(const DistTensor<T>& A, DistTensor<T>& B)
{
    DistTensor<T> *tmp = new DistTensor<T>(A);
  //  DistTensor<T> *tmp = new DistTensor<T>(A.TensorDist(), A.Grid());
  //  *tmp = A;

  for(Unsigned i = 0; i < A.Order(); ++i) {
    cout << i << " " << A.Dimension(i) << " vs. " << tmp->Dimension(i) << endl;
  }

  for (Unsigned mode = 0; mode < tmp->Order(); ++mode) {
    TensorDistribution dist = tmp->TensorDist();
    ModeDistribution modeDist = dist[mode];
    if (!(modeDist.empty())) {
      TensorDistribution newDist = dist;
      ModeDistribution newIgnoreDist = newDist[newDist.size() - 1];
      newIgnoreDist.insert(newIgnoreDist.end(), modeDist.begin(), modeDist.end());
      newDist[newDist.size() - 1] = newIgnoreDist;
      modeDist.clear();
      newDist[mode] = modeDist;
      DistTensor<T> *tmp2 = new DistTensor<T>(newDist, tmp->Grid());
      tmp2->GatherToOneRedistFrom(*tmp, mode);
      delete tmp;
      tmp = tmp2;
    }
  }

  if (TensorDistToString(B.TensorDist()) != TensorDistToString(tmp->TensorDist())) {
    cout << TensorDistToString(B.TensorDist()) << endl;
    cout << TensorDistToString(tmp->TensorDist()) << endl;
    throw;
  }
  B = *tmp;
  delete tmp;

  for(Unsigned i = 0; i < A.Order(); ++i) {
    cout << i << " " << A.Dimension(i) << " vs. " << B.Dimension(i) << endl;
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
	std::cout << "order " << i << " is " << gridDim << std::endl;
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

    ObjShape shapes3(3);
    shapes3[0] = 10;
    shapes3[1] = 10;
    shapes3[2] = 10;

    ObjShape shapes4(4);
    shapes4[0] = 10;
    shapes4[1] = 10;
    shapes4[2] = 10;
    shapes4[3] = 10;



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
    // B input has 4 dims
    //	Starting distribution: [D0,D1,D2,D3] or _D_0__D_1__D_2__D_3
    // C input has 3 dims
    //	Starting distribution: [D01,D2,D3] or _D_0_1__D_2__D_3


    DistTensor<T> A__D_0_1__D_2__D_3(shapes3, dist__D_0_1__D_2__D_3, g);
    DistTensor<T> B__D_0__D_1__D_2__D_3(shapes4, dist__D_0__D_1__D_2__D_3, g);
    DistTensor<T> C__D_0_1__D_2__D_3(shapes3, dist__D_0_1__D_2__D_3, g);

     Set(A__D_0_1__D_2__D_3);
     Set(B__D_0__D_1__D_2__D_3);
     Set(C__D_0_1__D_2__D_3);

     DistTensor<T> A_local( tmen::StringToTensorDist("[(),(),()]|(0,1,2,3)"), g );
     DistTensor<T> B_local( tmen::StringToTensorDist("[(),(),(),()]|(0,1,2,3)"), g );
     DistTensor<T> C_local( tmen::StringToTensorDist("[(),(),()]|(0,1,2,3)"), g );

     GatherAllModes(A__D_0_1__D_2__D_3, A_local);
     GatherAllModes(B__D_0__D_1__D_2__D_3, B_local);
     GatherAllModes(C__D_0_1__D_2__D_3, C_local);


    //**** (out of 4)
    //------------------------------------//

	tempShape = C__D_0_1__D_2__D_3.Shape();
	tempShape.push_back( g.Shape()[0] );
	tempShape.push_back( g.Shape()[3] );
	C__S__D_1__D_2__D_0__D_3.ResizeTo( tempShape );
	tempShape = C__D_0_1__D_2__D_3.Shape();
	C__D_1_0__D_2__D_3.ResizeTo( tempShape );
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
    C__D_1_0__D_2__D_3.ReduceScatterRedistFrom( C__D_1_0__D_2__S__D_3, 3, 2 );
    // C[D01,D2,D3] <- C[D10,D2,D3]
    C__D_0_1__D_2__D_3.PermutationRedistFrom( C__D_1_0__D_2__D_3, 0, modes_1_0 );

    //------------------------------------//

    //****


    IndexArray indices_aef( 3 );
    indices_aefcd[0] = 'a';
    indices_aefcd[1] = 'e';
    indices_aefcd[2] = 'f';

    for (int i = 0; i < 3; ++i) {
      cout << A_local.LockedTensor().Dimension(i) << endl;
    }
    cout << endl;
    for (int i = 0; i < 4; ++i) {
      cout << B_local.LockedTensor().Dimension(i) << endl;
    }
    cout << endl;
    for (int i = 0; i < 3; ++i) {
      cout << C_local.LockedTensor().Dimension(i) << endl;
    }
    cout << endl;

    LocalContract(1.0, A_local.LockedTensor(), indices_acd,
    		  B_local.LockedTensor(), indices_cefd,
    		  1.0, C_local.Tensor(), indices_aef);

    DistTensor<T> C_local_comparison( tmen::StringToTensorDist("[(),(),()]|(0,1,2,3)"), g );    
    GatherAllModes(C__D_0_1__D_2__D_3, C_local_comparison);

    /*
    Compare( C_local, C_local_comparison );
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
