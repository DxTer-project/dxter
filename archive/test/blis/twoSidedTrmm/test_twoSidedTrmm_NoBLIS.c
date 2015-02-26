#include <unistd.h>
#include "FLAME.h"

extern fla_axpy_t*  fla_axpy_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_hemm_t*  fla_hemm_cntl_blas;
extern fla_her2k_t* fla_her2k_cntl_blas;
extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
fla_eig_gest_t*     fla_eig_gest_ix_cntl;
fla_eig_gest_t*     fla_eig_gest_nx_cntl;
fla_eig_gest_t*     fla_eig_gest_ix_cntl_leaf;
fla_eig_gest_t*     fla_eig_gest_nx_cntl_leaf;

void FLA_TwoSidedTrmm( FLA_Obj A,
		       FLA_Obj B,
		       fla_eig_gest_t *cntl)
{
  FLA_Obj Y;

  dim_t m = FLA_Obj_length( A );
  

  FLA_Obj_create(FLA_DOUBLE,
		       m, m, 0, 0,
		       &Y);

  FLA_Eig_gest_nl( A, Y, B, cntl );
  //  err = FLA_Eig_gest_internal(FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, Y, B, 
  //			      fla_eig_gest_nx_cntl);
  //  err = FLA_Eig_gest(FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B);
  //err = FLA_Eig_gest(FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B);
  //    FLA_Eig_gest_nl_opt_var2( A, Y, B);

  FLA_Obj_free( &Y );
}


int main( int argc, char** argv )
{
  FLA_Obj a, b;
  FLA_Obj a_save, b_save;

  dim_t m, bs;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  dim_t bs_begin, bs_end, bs_inc;
  int   m_input, n_input, k_input;
  int   r, n_repeats;
  dim_t bsDxTBest, bsFLABest;

  double dtimeDxT;
  double dtime_saveDxT;
  double gflopsDxT;
  double gflopsDxTBest;
  double dtimeFLA;
  double dtime_saveFLA;
  double gflopsFLA;
  double gflopsFLABest;
  fla_eig_gest_t*     cntl;

  FLA_Init();

  n_repeats = 3;

#if 0
#define VAL 300;
  p_begin = VAL;
  p_end   = VAL;
  p_inc   = VAL;
  bs_begin = 256;
  bs_end   = 256;
  bs_inc = 16;
#else
  p_begin = 200;
  p_end   = 5000;
  p_inc   = 200;
  bs_begin = 64;
  bs_end   = 256;
  bs_inc   = 16;
#endif

  m_input = -1;


  for ( p = p_begin; p <= p_end; p += p_inc )
    {

      if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
      else               m =     ( dim_t )    m_input;

      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &a );
      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &a_save );
      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &b );
      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &b_save );

      FLA_Random_matrix( a );
      FLA_Random_matrix( b );

      FLA_Copy( a, a_save );
      FLA_Copy( b, b_save );

      gflopsDxTBest = 0;
      gflopsFLABest = 0;

      for ( bs = bs_begin; bs <= min(m,bs_end); bs += bs_inc ) {

	fla_blocksize_t *bs_fla = FLA_Blocksize_create(bs, bs, bs, bs);

	cntl = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
					     FLA_BLOCKED_VARIANT2,
					     bs_fla,
					     fla_eig_gest_nx_cntl_leaf,
					     fla_axpy_cntl_blas,
					     fla_axpy_cntl_blas,
					     fla_gemm_cntl_blas,
					     fla_gemm_cntl_blas,
					     fla_gemm_cntl_blas,
					     fla_hemm_cntl_blas,
					     fla_her2k_cntl_blas,
					     fla_trmm_cntl_blas,
					     fla_trmm_cntl_blas,
					     fla_trsm_cntl_blas,
					     fla_trsm_cntl_blas );

	dtime_saveDxT = 1.0e9;
	dtime_saveFLA = 1.0e9;

	for ( r = 0; r < n_repeats; ++r )
	  {
	    FLA_Copy( a_save, a );
	    FLA_Copy( b_save, b );

	    dtimeFLA = FLA_Clock();
	    FLA_TwoSidedTrmm( a, b, fla_eig_gest_nx_cntl);
	    dtime_saveFLA = min(dtimeFLA,FLA_Clock()-dtimeFLA);
	  }


	gflopsFLA = ( ((1.0 * m) * m) * m ) / ( dtime_saveFLA * 1.0e9 );
	if (gflopsFLA > gflopsFLABest) {
	  gflopsFLABest = gflopsFLA;
	  bsFLABest = bs;
	}
	
	printf( "data_TwoSidedTrmm_FLA_NoBLIS" );
	printf( "( %2lu, %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n\n",
		(p - p_begin + 1)/p_inc + 1, (bs - bs_begin + 1)/bs_inc + 1,
		m, bs, dtime_saveFLA, gflopsFLA );
	
	FLA_Cntl_obj_free( cntl );
	FLA_Blocksize_free( bs_fla );
      }

	printf( "data_TwoSidedTrmm_FLA_Best_NoBLIS" );
	printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
		(p - p_begin + 1)/p_inc + 1,
		m, gflopsFLABest, bsFLABest);
      printf("%%%%%%%%%%\n");
      fflush(stdout);

      FLA_Obj_free( &a );
      FLA_Obj_free( &b );
      FLA_Obj_free( &a_save );
      FLA_Obj_free( &b_save );
    }

  FLA_Finalize();

  return 0;
}

