#define INHERIT( par ) \
   INCLUDEPARENT( par )
#define __CLASSOP( pre, cla ) \
   if( type == pre##cla##_Type ) \
      return True;
#define _CLASSOP( pre, cla ) \
   __CLASSOP( pre, cla )
#undef CLASSOP
#define CLASSOP( cla ) \
   _CLASSOP( PREFIX, cla )
#define __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __VOIDABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames )
#define OVERRIDE( methodName, type, argTypes, argNames )
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )
#define __METHOD( pre, className, methodName, type, argTypes, argNames )
#define MEMBER( type, name )

#define _CLASSISSUPERBEGIN( pre, className ) \
   Bool pre##className##_IsSuper( Name type ) {
#define CLASSISSUPERBEGIN( pre, className ) \
   _CLASSISSUPERBEGIN( pre, className )
#define _CLASSISSUPEREND( pre, className ) \
      return False; \
   }
#define CLASSISSUPEREND( pre, className ) \
   _CLASSISSUPEREND( pre, className )

CLASSISSUPERBEGIN( PREFIX, CLASSNAME )
#include INCLUDECLASS( CLASSNAME )
CLASSISSUPEREND( PREFIX, CLASSNAME )

#undef INHERIT
#undef __CLASSOP
#undef _CLASSOP
#undef CLASSOP
#define CLASSOP( par )
#undef __VIRTUALMETHOD
#undef __VOIDVIRTUALMETHOD
#undef __ABSTRACTMETHOD
#undef __VOIDABSTRACTMETHOD
#undef OVERRIDE
#undef VOIDOVERRIDE
#undef __METHOD
#undef MEMBER
#undef _CLASSISSUPERBEGIN
#undef CLASSISSUPERBEGIN
#undef _CLASSISSUPEREND
#undef CLASSISSUPEREND
