#ifndef pcu_test_h
#define pcu_test_h

struct pcu_test_t {
      const char* name;
      pcu_suite_t* suite;
      pcu_testfunc_t* func;
      pcu_test_t* next;

      int globalresult;
      int nsrcs;
      pcu_source_t* srcs;
};

void pcu_test_run( pcu_test_t* test, pcu_listener_t* lsnr );
pcu_source_t* pcu_test_addSource( pcu_test_t* test, pcu_source_t* src );

#endif
