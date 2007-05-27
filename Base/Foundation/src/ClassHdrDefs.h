#define PARENT( classPath )				\
   INCLUDEFILE( StGermain/Base/Foundation, NoClass.h )
#define __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   typedef type (pre##className##_##methodName##Func)argTypes;
#define __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __VOIDABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define OVERRIDE( methodName, type, argTypes, argNames )
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )
#define __METHOD( pre, className, methodName, type, argTypes, argNames )
#define MEMBER( type, name )

#include INCLUDECLASS( CLASSNAME )

#undef PARENT
#undef __VIRTUALMETHOD
#undef __VOIDVIRTUALMETHOD
#undef __ABSTRACTMETHOD
#undef __VOIDABSTRACTMETHOD
#undef OVERRIDE
#undef VOIDOVERRIDE
#undef __METHOD
#undef MEMBER
