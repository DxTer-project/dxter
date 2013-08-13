#include <unistd.h>
#include "blis.h"
#include "SMBLAS3.h"

thread_comm_t global_comm[1];			
thread_comm_t proc_comms[NUMPROCS];		
thread_comm_t l2_comms[NUML2];
thread_comm_t l1_comms[NUML1];

void DxT_TrsmLLN( obj_t *alpha,
                 obj_t *L,
                 obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);
  FUNCTIONEND
}

void DxT_TrsmLLT( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);
  FUNCTIONEND
}


void DxT_TrsmLUN( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);
  FUNCTIONEND
}

void DxT_TrsmLUT( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);
  FUNCTIONEND
}

void DxT_TrsmRLN( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);
  FUNCTIONEND
}

void DxT_TrsmRLT( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);
  FUNCTIONEND
}

void DxT_TrsmRUN( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);
  FUNCTIONEND
}


void DxT_TrsmRUT( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);
  FUNCTIONEND
}

int main( int argc, char** argv )
{
  obj_t a, c1, c2;
  obj_t c_save;
  obj_t alpha, beta, negOne;
#ifdef TESTLIB
  obj_t normVal;
#endif
  dim_t m, n;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  int   m_input, n_input;
  num_t dt_a, dt_b, dt_c;
  num_t dt_alpha, dt_beta;
  int   r, n_repeats;

  double dtimeDxT;
  double dtime_saveDxT;
  double gflopsDxT;

  double dtimeFLA;
  double dtime_saveFLA;
  double gflopsFLA;

  int left=0, lower=0, trans=0;


  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;


  bli_obj_create( dt_alpha, 1, 1, 0, 0, &negOne );
  bli_setsc(  -1.0, 0.0, &negOne );		

  if (argc != 4) {
    printf("test L/R L/U N/T\n");
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

  if (*(argv[3]) == 'T')
    trans = 1;
  else if (*(argv[3]) != 'N') {
    printf("trans/normal not correct\n");
    return 0;
  }
  else
    trans = 0;



  bli_init();

  n_repeats = 3;

#if 1
  p_begin = 40;
  p_end   = 600;
  p_inc   = 40;
#else
  p_begin = 512;
  p_end   = 6144;
  p_inc   = 512;
#endif

  m_input = -1;
  n_input = -1;

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


		bli_obj_create( dt_alpha, 1, 1, 0, 0, &alpha );
		bli_obj_create( dt_beta,  1, 1, 0, 0, &beta );
#ifdef TESTLIB
		bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );
#endif
	
		if ( left )
			bli_obj_create( dt_a, m, m, 0, 0, &a );
		else
			bli_obj_create( dt_a, n, n, 0, 0, &a );
		bli_obj_create( dt_c, m, n, 0, 0, &c1 );
		bli_obj_create( dt_c, m, n, 0, 0, &c2 );
		bli_obj_create( dt_c, m, n, 0, 0, &c_save );


		bli_randm( &a );
		double *buff = (double*)bli_obj_buffer_at_off( a );
		dim_t rstride = bli_obj_row_stride( a );
		dim_t cstride = bli_obj_col_stride( a );
		for(unsigned int i=0;i<m;++i) {
		  buff[i*rstride+i*cstride] = 1;
		}
		bli_randm( &c1 );

		bli_setsc(  (1.0/1.0), 0.0, &alpha );
		bli_setsc(  (1.0/1.0), 0.0, &beta );
#ifdef TESTLIB
		bli_setsc(  (1.0/1.0), 0.0, &normVal );
#endif

		

		bli_copym( &c1, &c_save );
	
		dtime_saveDxT = 1.0e9;
		dtime_saveFLA = 1.0e9;

		for ( r = 0; r < n_repeats; ++r )
		{
			bli_copym( &c_save, &c1 );
			bli_copym( &c_save, &c2 );

			dtimeFLA = bli_clock();

			if (trans)
			  bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);

			bli_obj_set_struc( BLIS_TRIANGULAR, a );
			if (lower) {
			  bli_obj_set_uplo( BLIS_LOWER, a );
			}
			else {
			  bli_obj_set_uplo( BLIS_UPPER, a );
			}

			bli_trsm( left ? BLIS_LEFT : BLIS_RIGHT,
			          &alpha,
			          &a,
			          &c2 );

			if (trans)
			  bli_obj_set_conjtrans(BLIS_NO_TRANSPOSE, a);
			  // bli_obj_toggle_trans( a );
			bli_obj_set_struc( BLIS_GENERAL, a );
			bli_obj_set_uplo( BLIS_DENSE, a );
			
			dtime_saveFLA = bli_clock_min_diff( dtime_saveFLA, dtimeFLA );


			dtimeDxT = bli_clock();

			th_setup_comm(&global_comm[0], NUMTHREADS, 1);
			_Pragma( "omp parallel num_threads(NUMTHREADS)" ) 
			  {
			    if (left) {
			      if (lower) {
				if (trans) {
				  DxT_TrsmLLT(&alpha,
					      &a,
					      &c1 );
				}
				else { //normal
				  DxT_TrsmLLN(&alpha,
					      &a,
					      &c1 );
				}
			      }
			      else { //upper
				if (trans) {
				  DxT_TrsmLUT(&alpha,
					      &a,
					      &c1 );
				}
				else {//normal 
				  DxT_TrsmLUN(&alpha,
					      &a,
					      &c1 );
				}
			      }
			    }
			    else { //right
			      if (lower) {
				if (trans) {
				  DxT_TrsmRLT(&alpha,
					      &a,
					      &c1 );
				}
				else { //normal
				  DxT_TrsmRLN(&alpha,
					      &a,
					      &c1 );
				}

			      }
			      else { //upper
				if (trans) {
				  DxT_TrsmRUT(&alpha,
					      &a,
					      &c1 );
				}
				else { //normal
				  DxT_TrsmRUN(&alpha,
					      &a,
					      &c1 );
				}
			      }
			    }
			  }

			dtime_saveDxT = bli_clock_min_diff( dtime_saveDxT, dtimeDxT );
			/*
			bli_fnormm( &c1, &normVal );
			bli_printm( "c1 norm", &normVal, "%4.3e", "" );
			bli_fnormm( &c2, &normVal );
			bli_printm( "c2 norm", &normVal, "%4.3e", "" );
			*/

#ifdef TESTLIB
			bli_axpym( &negOne, &c1, &c2 );
			bli_fnormm( &c2, &normVal );
#endif

		}

		if ( left ) {
		  gflopsDxT = ( 1.0 * m * m * n ) / ( dtime_saveDxT * 1.0e9 );
		  gflopsFLA = ( 1.0 * m * m * n ) / ( dtime_saveFLA * 1.0e9 );
		}
		else {
		  gflopsDxT = ( 1.0 * m * n * n ) / ( dtime_saveDxT * 1.0e9 );
		  gflopsFLA = ( 1.0 * m * n * n ) / ( dtime_saveFLA * 1.0e9 );
		}

		printf( "data_trsm_DxT" );
		printf( "( %2ld, 1:4 ) = [ %4lu %4lu  %10.3e  %6.3f ];\n",
		        (p - p_begin + 1)/p_inc + 1, m, n, dtime_saveDxT, gflopsDxT );
		printf( "data_trsm_FLA" );
		printf( "( %2ld, 1:4 ) = [ %4lu %4lu  %10.3e  %6.3f ];\n\n",
		        (p - p_begin + 1)/p_inc + 1, m, n, dtime_saveFLA, gflopsFLA );

#ifdef TESTLIB
		bli_printm( "NORM", &normVal, "%4.1f", "" );
#endif
      

		bli_obj_free( &alpha );
		bli_obj_free( &beta );

		bli_obj_free( &a );
		bli_obj_free( &c1 );
		bli_obj_free( &c2 );
		bli_obj_free( &c_save );
	}

	bli_finalize();

	return 0;
}

