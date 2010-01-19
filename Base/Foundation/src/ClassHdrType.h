#define _CLASSTYPE( pre, className )		\
   extern const Name pre##className##_Type;
CLASSTYPE( PREFIX, CLASSNAME )

#undef _CLASSTYPE
