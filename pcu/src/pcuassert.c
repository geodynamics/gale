#include <stdlib.h>
#include "pcuassert.h"

int pcu_jump_ready = 0;
jmp_buf pcu_jump_env;
char* pcu_assert_cur = NULL;
