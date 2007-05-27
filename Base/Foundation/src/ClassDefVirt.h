#define PARENT( classPath )							\
  INCLUDEFILE( StGermain/Base/Foundation, NoClass.h )
#define __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
  type pre##className##_##methodName argTypes {					\
    assert( self );								\
    assert( ((pre##className*)self)->_##methodName );				\
    return ((pre##className*)self)->_##methodName argNames;			\
  }
#define __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
  type pre##className##_##methodName argTypes {					\
    assert( self );								\
    assert( ((pre##className*)self)->_##methodName );				\
    ((pre##className*)self)->_##methodName argNames;				\
  }
#define __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames )	\
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __VOIDABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define OVERRIDE( methodName, type, argTypes, argNames )	\
   _VIRTUALMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )	\
   _VOIDVIRTUALMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
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
