#define _CLASSTYPE( pre, className )			\
  const char* pre##className##_Type = #className;
CLASSTYPE( PREFIX, CLASSNAME )

#undef _CLASSTYPE
