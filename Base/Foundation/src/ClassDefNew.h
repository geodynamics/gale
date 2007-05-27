#define _CLASSNEW( pre, className )					\
  pre##className* pre##className##_New() {					\
    pre##className* self = (pre##className*)MemAlloc( pre##className, 		\
						      pre##className##_Type );	\
    pre##className##_Construct( self );						\
    return self;								\
  }
#define CLASSNEW( pre, className )		\
   _CLASSNEW( pre, className )

CLASSNEW( PREFIX, CLASSNAME )

#undef _CLASSNEW
#undef CLASSNEW
