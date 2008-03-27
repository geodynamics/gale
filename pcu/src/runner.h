#ifndef pcu_runner_h
#define pcu_runner_h

void pcu_runner_init( int argc, char* argv[] );
void pcu_runner_finalise();
void pcu_runner_run( pcu_listener_t* lsnr );
void _pcu_runner_addSuite( const char* name,
			   void (initfunc)( pcu_suite_t* ) );

#define pcu_runner_addSuite( suite )            \
   _pcu_runner_addSuite( #suite, suite )

#endif
