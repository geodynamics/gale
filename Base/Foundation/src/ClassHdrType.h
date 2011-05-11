#define _CLASSTYPE( pre, className )		\
   extern const Type pre##className##_Type;
CLASSTYPE( PREFIX, CLASSNAME )

#undef _CLASSTYPE
