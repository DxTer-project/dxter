#include <unistd.h>
#include "FLAME.h"

fla_qrut_t *qrut_cntl;

//#define TESTLIB
void FLA_QR( FLA_Obj A, 
	     FLA_Obj T )
{
  FLA_QR_UT( A, T );
}


int main( int argc, char** argv )
{
  FLA_Obj a;
  FLA_Obj a_save;
  FLA_Obj t;
  FLA_Obj t_save;
#ifdef TESTLIB
  obj_t normValA;
#endif
  dim_t m;
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

  FLA_Init();

  qrut_cntl = FLA_Cntl_qrut_obj_create( FLA_FLAT,
					FLA_UNB_OPT_VARIANT2,
					NULL,
					NULL,
					NULL ); 

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
      FLA_Obj_create( FLA_DOUBLE, m, m, 0, 0, &a_save );

      FLA_Random_matrix( a );

      FLA_Copy( a, a_save );

      gflopsDxTBest = 0;
      gflopsFLABest = 0;

      for ( dim_t bs = bs_begin; bs <= bs_end; bs += bs_inc ) {
	FLA_Obj_create( FLA_DOUBLE, bs, m, 0, 0, &t );
	FLA_Obj_create( FLA_DOUBLE, bs, m, 0, 0, &t_save );
	FLA_Random_matrix( t );
	FLA_Copy( t, t_save );
	
	
	dtime_saveDxT = 1.0e9;
	dtime_saveFLA = 1.0e9;

	for ( r = 0; r < n_repeats; ++r )
	  {
	    FLA_Copy( a_save, a );
	    FLA_Copy( t_save, t );

	    dtimeFLA = FLA_Clock();
	    FLA_QR( a, t );
	    dtime_saveFLA = min(FLA_Clock()-dtimeFLA,dtimeFLA);
	  }

	gflopsFLA = ( (4.0 / 3.0) * m * m * m ) / ( dtime_saveFLA * 1.0e9 );
	if (gflopsFLA > gflopsFLABest) {
	  gflopsFLABest = gflopsFLA;
	  bsFLABest = bs;
	}

	printf( "data_QR_FLA_NoBLIS" );
	printf( "( %2lu, %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n\n",
		(p - p_begin + 1)/p_inc + 1, (bs - bs_begin + 1)/bs_inc + 1,
		m, bs, dtime_saveFLA, gflopsFLA );

	FLA_Obj_free( &t );
	FLA_Obj_free( &t_save );
      }
      
      printf( "data_QR_FLA_Best_NoBLIS" );
      printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
	      (p - p_begin + 1)/p_inc + 1,
	      m, gflopsFLABest, bsFLABest);

      printf("%%%%%%%%%%\n");
      fflush(stdout);

      FLA_Obj_free( &a );
      FLA_Obj_free( &a_save );
    }

  FLA_Finalize();

  return 0;
}

