#include <unistd.h>
#include "FLAME.h"

extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_appiv_t* fla_appiv_cntl_leaf;

fla_lu_t*           fla_lu_piv_cntl;
fla_lu_t*           fla_lu_piv_cntl_in;
fla_lu_t*           fla_lu_piv_cntl_leaf;
fla_blocksize_t*    fla_lu_piv_var5_bsize;
fla_blocksize_t*    fla_lu_piv_var5_bsize_in;

void FLA_LU( FLA_Obj A, 
	     FLA_Obj p,
	     fla_lu_t* fla_lu_piv_cntl2)
{
  FLA_LU_piv_internal( A, p, fla_lu_piv_cntl2 );
}


int main( int argc, char** argv )
{
  FLA_Obj pVec, a;
  FLA_Obj a_save;
  dim_t m, bs;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  dim_t bs_begin, bs_end, bs_inc;
  int   m_input, n_input, k_input;
  int   r, n_repeats;

  double dtimeDxT;
  double dtime_saveDxT;
  double gflopsDxT;
  double gflopsDxTBest;
  double dtimeFLA;
  double dtime_saveFLA;
  double gflopsFLA;
  double gflopsFLABest;
  dim_t bsDxTBest, bsFLABest;

  fla_lu_t*           fla_lu_piv_cntl2;


  FLA_Init();

  n_repeats = 3;

  p_begin = 200;
  p_end   = 5000;
  p_inc   = 200;
  bs_begin = 64;
  bs_end   = 256;
  bs_inc   = 16;

  m_input = -1;

  for ( p = p_begin; p <= p_end; p += p_inc )
    {

      if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
      else               m =     ( dim_t )    m_input;

      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &a );
      FLA_Obj_create( FLA_INT, m, 1, 0, 0, &pVec );
      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &a_save );

      FLA_Random_matrix( a );

      FLA_Copy( a, a_save );
	
      gflopsDxTBest = 0;
      gflopsFLABest = 0;

      for ( bs = bs_begin; bs <= bs_end; bs += bs_inc ) {

	fla_blocksize_t *bs_fla = FLA_Blocksize_create(bs, bs, bs, bs);


	fla_lu_piv_cntl2       = FLA_Cntl_lu_obj_create( FLA_FLAT,
                                                         FLA_BLOCKED_VARIANT5,
							 bs_fla,
                                                         fla_lu_piv_cntl_in,
                                                         fla_gemm_cntl_blas,
                                                         fla_gemm_cntl_blas,
                                                         fla_gemm_cntl_blas,
                                                         fla_trsm_cntl_blas,
                                                         fla_trsm_cntl_blas,
                                                         fla_appiv_cntl_leaf,
                                                         fla_appiv_cntl_leaf );




	dtime_saveDxT = 1.0e9;
	dtime_saveFLA = 1.0e9;

	for ( r = 0; r < n_repeats; ++r )
	  {
	    FLA_Copy( a_save, a );

	    dtimeFLA = FLA_Clock();
	    FLA_LU( a, pVec, fla_lu_piv_cntl2 );
	    dtime_saveFLA = min(dtimeFLA,FLA_Clock()-dtimeFLA);
	  }

	gflopsFLA = ( (2.0 / 3.0) * m * m * m ) / ( dtime_saveFLA * 1.0e9 );
	if (gflopsFLA > gflopsFLABest) {
	  gflopsFLABest = gflopsFLA;
	  bsFLABest = bs;
	}

	printf( "data_LU_FLA_NoBLIS" );
	printf( "( %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n",
		(p - p_begin + 1)/p_inc + 1, m, bs, dtime_saveFLA, gflopsFLA );

	FLA_Cntl_obj_free( fla_lu_piv_cntl2 );
	FLA_Blocksize_free( bs_fla );
      }

	printf( "data_LU_FLA_Best_NoBLIS" );
	printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
		(p - p_begin + 1)/p_inc + 1,
		m, gflopsFLABest, bsFLABest);
      printf("%%%%%%%%%%\n\n\n");



      FLA_Obj_free( &a );
      FLA_Obj_free( &pVec );
      FLA_Obj_free( &a_save );
    }

  FLA_Finalize();

  return 0;
}

