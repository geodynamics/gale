#define INHERIT( par ) \
   INCLUDEPARENT( par )
#define __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   pre##className##_##methodName##Func* _##methodName;
#define __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __VOIDABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define OVERRIDE( methodName, type, argTypes, argNames )
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )
#define __METHOD( pre, className, methodName, type, argTypes, argNames )
#define MEMBER( type, name ) \
   type	name;
#define _CLASSDECL( pre, name ) \
   struct pre##name
#define CLASSDECL( pre, name ) \
   _CLASSDECL( pre, name )

CLASSDECL( PREFIX, CLASSNAME ) {
#include INCLUDECLASS( CLASSNAME )
};

#undef INHERIT
#undef __VIRTUALMETHOD
#undef __VOIDVIRTUALMETHOD
#undef __ABSTRACTMETHOD
#undef __VOIDABSTRACTMETHOD
#undef OVERRIDE
#undef VOIDOVERRIDE
#undef __METHOD
#undef MEMBER
#undef CLASSDECL
