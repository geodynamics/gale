#ifndef pcu_listener_h
#define pcu_listener_h

struct pcu_listener_t {
      pcu_suiteentry_t* suitebegin;
      pcu_suiteentry_t* suiteend;
      pcu_testentry_t* testbegin;
      pcu_testentry_t* testend;
      pcu_assertentry_t* assertdone;
      void* data;
};

#endif
