#include <unistd.h>
#include "FLAME.h"


extern fla_axpy_t*  fla_axpy_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_hemm_t*  fla_hemm_cntl_blas;
extern fla_her2k_t* fla_her2k_cntl_blas;
extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;

fla_eig_gest_t*     fla_eig_gest_ix_cntl_leaf;
fla_eig_gest_t*     fla_eig_gest_nx_cntl_leaf;

void FLA_TwoSidedTrsm( FLA_Obj A, 
		       FLA_Obj B,
		       fla_eig_gest_t *cntl )
{
  FLA_Obj Y;

  FLA_Obj_create( FLA_DOUBLE, 	 
		  FLA_Obj_length(A), FLA_Obj_width(B),
		  0, 0, &Y );
  
#if 1
  FLA_Eig_gest_il( A, Y, B, cntl );
#elif 1
  FLA_Eig_gest(FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B);
#else
  FLA_Eig_gest_il_opt_var2( A, Y, B);
#endif

  FLA_Obj_free(&Y);
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
  double gflopsDxT, gflopsDxTBest;
  double dtimeFLA;
  double dtime_saveFLA;
  double gflopsFLA, gflopsFLABest;


  fla_eig_gest_t*     cntl;

  FLA_Init();


  n_repeats = 1;

#if 0
#define VAL 130
  p_begin = VAL;
  p_end   = VAL;
  p_inc   = VAL;
#else
  p_begin = 200;
  p_end   = 5000;
  p_inc   = 200;
  bs_begin = 64;
  bs_end   = 256;
  bs_inc   = 16;
#endif

  m_input = -1;

      //  bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );

  for ( p = p_begin; p <= p_end; p += p_inc )
    {

      if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
      else               m =     ( dim_t )    m_input;

      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &a );
      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &a_save );
      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &b );
      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &b_save );
      
      FLA_Random_matrix( a );
      FLA_Random_herm_matrix( FLA_LOWER_TRIANGULAR, b );

      FLA_Copy( a, a_save );
      FLA_Copy( b, b_save );
	
      gflopsDxTBest = 0;
      gflopsFLABest = 0;
      
      for ( bs = bs_begin; bs <= bs_end; bs += bs_inc ) {

	fla_blocksize_t *bs_fla = FLA_Blocksize_create(bs, bs, bs, bs);

	cntl = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
					     FLA_BLOCKED_VARIANT4,
					     bs_fla,
					     fla_eig_gest_ix_cntl_leaf,
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
	    FLA_TwoSidedTrsm( a, b, cntl);
	    dtime_saveFLA = min(dtimeFLA, FLA_Clock()-dtimeFLA);
	  }

	gflopsFLA = ( ((1.0 * m) * m) * m ) / ( dtime_saveFLA * 1.0e9 );
	if (gflopsFLA > gflopsFLABest) {
	  gflopsFLABest = gflopsFLA;
	  bsFLABest = bs;
	}
	
	printf( "data_TwoSidedTrsm_FLA_NoBLIS" );
	printf( "( %2lu, %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n\n",
		(p - p_begin + 1)/p_inc + 1, (bs - bs_begin + 1)/bs_inc + 1,
		m, bs, dtime_saveFLA, gflopsFLA );

	fflush(stdout);
	
	FLA_Cntl_obj_free( cntl );
	FLA_Blocksize_free( bs_fla );
      }

	printf( "data_TwoSidedTrsm_FLA_Best_NoBLIS" );
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

