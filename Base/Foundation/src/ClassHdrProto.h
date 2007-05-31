#define _CLASSNEW( pre, className ) \
   pre##className* pre##className##_New();
#define CLASSNEW( pre, className ) \
   _CLASSNEW( pre, className )
#define _CLASSCONSTRUCT( pre, className ) \
   void pre##className##_Construct( void* self );
#define CLASSCONSTRUCT( pre, className ) \
   _CLASSCONSTRUCT( pre, className )
#define _CLASSISSUPER( pre, className ) \
   Bool pre##className##_IsSuper( const char* type );
#define CLASSISSUPER( pre, className ) \
   _CLASSISSUPER( pre, className )

CLASSNEW( PREFIX, CLASSNAME )
CLASSCONSTRUCT( PREFIX, CLASSNAME )
CLASSISSUPER( PREFIX, CLASSNAME )

#undef _CLASSNEW
#undef CLASSNEW
#undef _CLASSCONSTRUCT
#undef CLASSCONSTRUCT
#undef _CLASSISSUPER
#undef CLASSISSUPER

#define INHERIT( par ) \
   INCLUDEFILE( StGermain/Base/Foundation, ClassEmpty.h )
#define __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   type _##pre##className##_##methodName argTypes;				\
   type pre##className##_##methodName argTypes;
#define __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   type pre##className##_##methodName argTypes;
#define __VOIDABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames )
#define OVERRIDE( methodName, type, argTypes, argNames )	\
   _VIRTUALMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )	\
   _VOIDVIRTUALMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define __METHOD( pre, className, methodName, type, argTypes, argNames )	\
   type pre##className##_##methodName argTypes;
#define MEMBER( type, name )

#include INCLUDECLASS( CLASSNAME )

#undef INHERIT
#undef __VIRTUALMETHOD
#undef __VOIDVIRTUALMETHOD
#undef __ABSTRACTMETHOD
#undef __VOIDABSTRACTMETHOD
#undef OVERRIDE
#undef VOIDOVERRIDE
#undef __METHOD
#undef MEMBER
