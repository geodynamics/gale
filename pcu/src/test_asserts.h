#ifndef test_asserts_h
#define test_asserts_h

#include "types.h"
#include "listener.h"
#include "suite.h"
#include "source.h"

extern pcu_suite_t* pcu_cursuite;

#define _pcu_assert_eval( expr, type )                  \
   pcu_cursuite->lsnr->assertdone(                      \
      pcu_cursuite->lsnr,                               \
      pcu_test_addSource(                               \
	 pcu_cursuite->curtest,                         \
	 pcu_source_create( (expr) ? 1 : 0,             \
                            type,			\
                            __FILE__,                   \
                            __LINE__,                   \
                            #expr,			\
                            NULL,			\
                            pcu_cursuite->curtest )     \
	 )                                              \
      )

#define pcu_assert_true( expr )			\
   _pcu_assert_eval( expr, "true" )

#define pcu_assert_gt( a, b )				\
   _pcu_assert_eval( (a) > (b), "greater than" )

#define pcu_assert_lt( a, b )			\
   _pcu_assert_eval( (a) < (b), "less than" )

#define pcu_assert_ge( a, b )					\
   _pcu_assert_eval( (a) >= (b), "greater than or equal" )

#define pcu_assert_le( a, b )				\
   _pcu_assert_eval( (a) <= (b), "less than or equal" )

#endif
