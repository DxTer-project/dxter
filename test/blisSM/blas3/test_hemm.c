#include <unistd.h>
#include "blis.h"
#include "SMBLAS3.h"

thread_comm_t global_comm[1];			
thread_comm_t proc_comms[NUMPROCS];		
thread_comm_t l2_comms[NUML2];
thread_comm_t l1_comms[NUML1];


void DxT_HemmLL( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  FUNCTIONSTART
  bli_scalm(beta, C);
  FUNCTIONEND
}

void DxT_HemmLU( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  FUNCTIONSTART
  bli_scalm(beta, C);
  FUNCTIONEND

}

void DxT_HemmRL( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  FUNCTIONSTART
  bli_scalm(beta, C);
  FUNCTIONEND

}

void DxT_HemmRU( obj_t *alpha,
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
	obj_t alpha, beta, normVal;
	dim_t m, n;
	dim_t p;
	dim_t p_begin, p_end, p_inc;
	int   m_input, n_input;
	num_t dt_a, dt_b, dt_c;
	num_t dt_alpha, dt_beta;
	int   r, n_repeats;
	side_t side;

	double dtime;
	double dtime_save;
	double gflops;

	int left, lower;

	bli_init();

	n_repeats = 3;

	p_begin = 40;
	p_end   = 600;
	p_inc   = 40;

	m_input = -1;
	n_input = -1;

	dt_a = BLIS_DOUBLE;
	dt_b = BLIS_DOUBLE;
	dt_c = BLIS_DOUBLE;
	dt_alpha = BLIS_DOUBLE;
	dt_beta = BLIS_DOUBLE;

	if (argc != 3) {
	  printf("test L/R L/U\n");
	  fflush(stdout);
	  return 0;
	}

	if (*(argv[1]) == 'L')
	  left = 1;
	else if (*(argv[1]) != 'R') {
	  printf("left/right not correct\n");
	  return 0;
	}
	else
	  left = 0;

	if (*(argv[2]) == 'L')
	  lower = 1;
	else if (*(argv[2]) != 'U') {
	  printf("lower/upper not correct\n");
	  return 0;
	}
	else
	  lower = 0;


	for ( p = p_begin; p <= p_end; p += p_inc )
	{

		if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
		else               m =     ( dim_t )    m_input;
		if ( n_input < 0 ) n = p * ( dim_t )abs(n_input);
		else               n =     ( dim_t )    n_input;


		bli_obj_create( dt_alpha, 1, 1, 0, 0, &alpha );
		bli_obj_create( dt_beta,  1, 1, 0, 0, &beta );
		bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );
	

		if ( left )
		  bli_obj_create( dt_a, m, m, 0, 0, &a );
		else
		  bli_obj_create( dt_a, n, n, 0, 0, &a );
		bli_obj_create( dt_b, m, n, 0, 0, &b );
		bli_obj_create( dt_c, m, n, 0, 0, &c1 );
		bli_obj_create( dt_c, m, n, 0, 0, &c2 );
		bli_obj_create( dt_c, m, n, 0, 0, &c_save );

		bli_randm( &a );
		bli_randm( &b );
		bli_randm( &c1 );


		bli_setsc( 1.0, 0.0, &alpha );
		bli_setsc( 1.0, 0.0, &beta );

		bli_copym( &c1, &c_save );
	
		dtime_save = 1.0e9;

		for ( r = 0; r < n_repeats; ++r )
		{
			bli_copym( &c_save, &c1 );
			bli_copym( &c_save, &c2 );


			dtime = bli_clock();

			bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );
			
			bli_obj_set_struc( BLIS_HERMITIAN, a );
			if (lower) {
			  bli_obj_set_uplo( BLIS_LOWER, a );
			}
			else {
			  bli_obj_set_uplo( BLIS_UPPER, a );
			}
			
			bli_hemm( left ? BLIS_LEFT : BLIS_RIGHT,
				  &alpha,
				  &a,
				  &b,
				  &beta,
				  &c1);

			bli_obj_set_struc( BLIS_GENERAL, a );
			bli_obj_set_uplo( BLIS_DENSE, a );


			th_setup_comm(&global_comm[0], NUMTHREADS, 1);
			_Pragma( "omp parallel num_threads(NUMTHREADS)" ) 
			  {
			    if (left) {
			      if (lower) {
				DxT_HemmLL( &alpha,
					    &a,
					    &b,
					    &beta,
					    &c2);
			      }
			      else { //upper
				DxT_HemmLU( &alpha,
					    &a,
					    &b,
					    &beta,
					    &c2);
			      }
			    }
			    else { //right
			      if (lower) {
				DxT_HemmRL( &alpha,
					    &a,
					    &b,
					    &beta,
					    &c2);
			      }
			      else { //upper
				DxT_HemmRU( &alpha,
					    &a,
					    &b,
					    &beta,
					    &c2);
			      }
			    }
			  }

			bli_axpym( &BLIS_MINUS_ONE, &c1, &c2 );
			
			bli_fnormm( &c2, &normVal );

			dtime_save = bli_clock_min_diff( dtime_save, dtime );
		}

		if ( left )
			gflops = ( 2.0 * m * m * n ) / ( dtime_save * 1.0e9 );
		else
			gflops = ( 2.0 * m * n * n ) / ( dtime_save * 1.0e9 );

		printf( "data_hemm_blis" );
		printf( "( %2ld, 1:4 ) = [ %4lu %4lu  %10.3e  %6.3f ];\n",
		        (p - p_begin + 1)/p_inc + 1, m, n, dtime_save, gflops );

		bli_printm( "NORM", &normVal, "%4.1f", "" );

		bli_obj_free( &alpha );
		bli_obj_free( &beta );
		bli_obj_free( &normVal );

		bli_obj_free( &a );
		bli_obj_free( &b );
		bli_obj_free( &c1 );
		bli_obj_free( &c2 );
		bli_obj_free( &c_save );
	}

	bli_finalize();

	return 0;
}

