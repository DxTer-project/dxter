#include <unistd.h>
#include <omp.h>
#include "blis.h"
#include "SMBLAS3.h"

thread_comm_t global_comm[1];
thread_comm_t all_l2_comms[1];
thread_comm_t proc_comms[NUMPROCS];
thread_comm_t l2_comms[NUML2];
thread_comm_t l2_comms_sub[NUML2];
thread_comm_t l1_comms[NUML1];


void DxT_GemmNN( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  FUNCTIONSTART

  bli_scalm(beta, C);



  FUNCTIONEND
}

void DxT_GemmTN( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  FUNCTIONSTART

  bli_scalm(beta, C);
  
  FUNCTIONEND
}

void DxT_GemmNT( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  FUNCTIONSTART

  bli_scalm(beta, C);

  FUNCTIONEND
}

void DxT_GemmTT( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  FUNCTIONSTART

  bli_scalm(beta, C);

FUNCTIONEND
}

int main( int argc, char** argv )
{
  obj_t a, b, c1, c2;
  obj_t c_save;
  obj_t alpha, beta, negOne, normVal;
  dim_t m, n, k;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  int   m_input, n_input, k_input;
  num_t dt_a, dt_b, dt_c;
  num_t dt_alpha, dt_beta;
  int   r, n_repeats;

  double dtime;
  double dtime_save;
  double gflops;

  int transA=0, transB=0;

  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;

  bli_obj_create( dt_alpha, 1, 1, 0, 0, &negOne );
  bli_setsc(  -1.0, 0.0, &negOne );		
			
  if (argc != 3) {
    printf("test T/N T/N\n");
    return 0;
  }

  if (*(argv[1]) == 'T')
    transA = 1;
  else if (*(argv[1]) != 'N') {
    printf("transA not correct\n");
    return 0;
  }

  if (*(argv[2]) == 'T')
    transB = 1;
  else if (*(argv[2]) != 'N') {
    printf("transB not correct\n");
    return 0;
  }

  bli_init();

  n_repeats = 1;

  p_begin = 4;
  p_end   = 12;
  p_inc   = 4;

  m_input = -1;
  //m_input = 384;
  n_input = -1;
  k_input = -1;
  //k_input = 200;

  dt_a = BLIS_DOUBLE;
  dt_b = BLIS_DOUBLE;
  dt_c = BLIS_DOUBLE;
  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;

  for ( p = p_begin; p <= p_end; p += p_inc )
    {

      if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
      else               m =     ( dim_t )    m_input;
      if ( n_input < 0 ) n = p * ( dim_t )abs(n_input);
      else               n =     ( dim_t )    n_input;
      if ( k_input < 0 ) k = p * ( dim_t )abs(k_input);
      else               k =     ( dim_t )    k_input;


      bli_obj_create( dt_alpha, 1, 1, 0, 0, &alpha );
      bli_obj_create( dt_beta,  1, 1, 0, 0, &beta );

      bli_obj_create( dt_a, m, k, 0, 0, &a );
      bli_obj_create( dt_b, k, n, 0, 0, &b );
      bli_obj_create( dt_c, m, n, 0, 0, &c1 );
      bli_obj_create( dt_c, m, n, 0, 0, &c2 );
      bli_obj_create( dt_c, m, n, 0, 0, &c_save );

      bli_randm( &a );
      bli_randm( &b );
      bli_randm( &c1 );


      bli_setsc( 1.0, 0.0, &alpha );
      bli_setsc( 1.0/1.0, 0.0, &beta );

      bli_copym( &c1, &c_save );
	
      dtime_save = 1.0e9;

      for ( r = 0; r < n_repeats; ++r )
	{
	  bli_copym( &c_save, &c1 );
	  bli_copym( &c_save, &c2 );


	  dtime = bli_clock();

	  //bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );
	  if (!transA && !transB) {
	    th_setup_comm(&global_comm[0], NUMTHREADS, 1);
	    //	    _Pragma( "omp parallel num_threads(gemm_num_threads_default)" ) 
	    _Pragma( "omp parallel num_threads(NUMTHREADS)" ) 
	      {
		DxT_GemmNN( &alpha,
			    &a,
			    &b,
			    &beta,
			    &c2 );
	      }
			
	    bli_gemm( &alpha,
		      &a,
		      &b,
		      &beta,
		      &c1);
	  }
	  else if (transA && !transB) {
	    DxT_GemmTN( &alpha,
			&a,
			&b,
			&beta,
			&c2 );

	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);
	    bli_gemm( &alpha,
		      &a,
		      &b,
		      &beta,
		      &c1);
	    bli_obj_toggle_trans( a );
	  }
	  else if (!transA) {// &&transB 
	    DxT_GemmNT( &alpha,
			&a,
			&b,
			&beta,
			&c2 );

	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, b);
	    bli_gemm( &alpha,
		      &a,
		      &b,
		      &beta,
		      &c1);
	    bli_obj_toggle_trans( b );
	  }
	  else {//transA && transB
	    DxT_GemmTT( &alpha,
			&a,
			&b,
			&beta,
			&c2 );

	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);
	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, b);
	    bli_gemm( &alpha,
		      &a,
		      &b,
		      &beta,
		      &c1);
	    bli_obj_toggle_trans( a );
	    bli_obj_toggle_trans( b );
	  }
	  
	  bli_printm( "c_save", &c_save, "%4.1f", "" );
	  bli_printm( "c1", &c1, "%4.1f", "" );
	  bli_printm( "c2", &c2, "%4.1f", "" );

	  bli_axpym( &negOne, &c1, &c2 );

	  bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );

	  bli_fnormm( &c2, &normVal );

	  dtime_save = bli_clock_min_diff( dtime_save, dtime );
	}

      gflops = ( 2.0 * m * k * n ) / ( dtime_save * 1.0e9 );

      printf( "data_gemm_blis" );
      printf( "( %2ld, 1:5 ) = [ %4lu %4lu %4lu  %10.3e  %6.3f ];\n",
	      (p - p_begin + 1)/p_inc + 1, m, k, n, dtime_save, gflops );

      bli_printm( "NORM", &normVal, "%4.1f", "" );

      bli_obj_free( &alpha );
      bli_obj_free( &beta );

      bli_obj_free( &a );
      bli_obj_free( &b );
      bli_obj_free( &c1 );
      bli_obj_free( &c2 );
      bli_obj_free( &c_save );
    }

  bli_finalize();

  return 0;
}

