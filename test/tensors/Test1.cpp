/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "tensormental.hpp" instead
#include "tensormental.hpp"
using namespace tmen;

#define GRIDORDER 4

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
    ModeArray modes_1;
    modes_1.push_back(1);
    ModeArray modes_2;
    modes_2.push_back(2);
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

    DistTensor<T> A(shapes3, dist__D_0_1__D_2__D_3, indices_acd, g);
    DistTensor<T> B(shapes4, dist__D_0__D_1__D_2__D_3, indices_cefd, g);
    DistTensor<T> C(shapes3, dist__D_0_1__D_2__D_3, indices_acd, g);

    Set(A);
    Set(B);
    Set(C);

    //**** (out of 4)
    //------------------------------------//

    //**** (out of 2)
    //------------------------------------//

    // A[D0,D2,D3] <- A[D01,D2,D3]
    A__D_0__D_2__D_3.ResizeTo( A__D_0_1__D_2__D_3 );
    AllGatherRedist( A__D_0__D_2__D_3, A__D_0_1__D_2__D_3, 0, modes_1 );
    // A[*,D2,D3] <- A[D0,D2,D3]
    A__S__D_2__D_3.ResizeTo( A__D_0__D_2__D_3 );
    AllGatherRedist( A__S__D_2__D_3, A__D_0__D_2__D_3, 0, modes_0 );
    // A[*,D20,D3] <- A[*,D2,D3]
    A__S__D_2_0__D_3.ResizeTo( A__S__D_2__D_3 );
    LocalRedist( A__S__D_2_0__D_3, A__S__D_2__D_3, 1, modes_0 );
    // A[*,D02,D3] <- A[*,D20,D3]
    A__S__D_0_2__D_3.ResizeTo( A__S__D_2_0__D_3 );
    PermutationRedist( A__S__D_0_2__D_3, A__S__D_2_0__D_3, 1 );
    // A[*,D0,D3] <- A[*,D02,D3]
    A__S__D_0__D_3.ResizeTo( A__S__D_0_2__D_3 );
    AllGatherRedist( A__S__D_0__D_3, A__S__D_0_2__D_3, 1, modes_2 );

    //------------------------------------//

    //****
    // 1.0 * A[*,D0,D3]_acd * B[D0,D1,D2,D3]_cefd + 0.0 * C[*,D1,D2,D0,D3]_aefcd
    A__S__D_0__D_3.SetIndices( indices_acd );
    B__D_0__D_1__D_2__D_3.SetIndices( indices_cefd );
    C__S__D_1__D_2__D_0__D_3.SetIndices( indices_aefcd );
    LocalContract(1.0, A__S__D_0__D_3.LockedTensor(), B__D_0__D_1__D_2__D_3.LockedTensor(), 0.0, C__S__D_1__D_2__D_0__D_3.Tensor());
    // C[*,D1,*,D0,D3] <- C[*,D1,D2,D0,D3]
    C__S__D_1__S__D_0__D_3.ResizeTo( C__S__D_1__D_2__D_0__D_3 );
    AllGatherRedist( C__S__D_1__S__D_0__D_3, C__S__D_1__D_2__D_0__D_3, 2, modes_2 );
    // C[*,D12,*,D0,D3] <- C[*,D1,*,D0,D3]
    C__S__D_1_2__S__D_0__D_3.ResizeTo( C__S__D_1__S__D_0__D_3 );
    LocalRedist( C__S__D_1_2__S__D_0__D_3, C__S__D_1__S__D_0__D_3, 1, modes_2 );
    // C[*,D21,*,D0,D3] <- C[*,D12,*,D0,D3]
    C__S__D_2_1__S__D_0__D_3.ResizeTo( C__S__D_1_2__S__D_0__D_3 );
    PermutationRedist( C__S__D_2_1__S__D_0__D_3, C__S__D_1_2__S__D_0__D_3, 1 );
    // C[*,D2,*,D0,D3] <- C[*,D21,*,D0,D3]
    C__S__D_2__S__D_0__D_3.ResizeTo( C__S__D_2_1__S__D_0__D_3 );
    AllGatherRedist( C__S__D_2__S__D_0__D_3, C__S__D_2_1__S__D_0__D_3, 1, modes_1 );
    // C[D1,D2,*,D0,D3] <- C[*,D2,*,D0,D3]
    C__D_1__D_2__S__D_0__D_3.ResizeTo( C__S__D_2__S__D_0__D_3 );
    LocalRedist( C__D_1__D_2__S__D_0__D_3, C__S__D_2__S__D_0__D_3, 0, modes_1 );
    // C[D10,D2,*,D3] <- C[D1,D2,*,D0,D3] (with SumScatter on D0)
    C__D_1_0__D_2__S__D_3.ResizeTo( C__D_1__D_2__S__D_0__D_3 );
    ReduceScatterRedist( C__D_1_0__D_2__S__D_3, C__D_1__D_2__S__D_0__D_3, 3, 0 );
    // C[D10,D2,D3] <- C[D10,D2,*,D3] (with SumScatter on D3)
    C__D_1_0__D_2__D_3.ResizeTo( C__D_1_0__D_2__S__D_3 );
    ReduceScatterRedist( C__D_1_0__D_2__D_3, C__D_1_0__D_2__S__D_3, 3, 2 );
    // C[D01,D2,D3] <- C[D10,D2,D3]
    C__D_0_1__D_2__D_3.ResizeTo( C__D_1_0__D_2__D_3 );
    PermutationRedist( C__D_0_1__D_2__D_3, C__D_1_0__D_2__D_3, 0 );

    //------------------------------------//

    //****


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

        const Grid g( comm, GRIDORDER, args.gridShape );

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
