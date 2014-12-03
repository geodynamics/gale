#define _CLASSTYPE( pre, className )			\
  const Type pre##className##_Type = #className;
CLASSTYPE( PREFIX, CLASSNAME )

#undef _CLASSTYPE
