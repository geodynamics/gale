#define _CLASSTYPE( pre, className )		\
   extern const char* pre##className##_Type;
CLASSTYPE( PREFIX, CLASSNAME )

#undef _CLASSTYPE
