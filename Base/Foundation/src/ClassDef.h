#ifndef CLASSNAME
#error No class name given.
#endif

#include "ClassPre.h"

#define _CLASSTYPE( className )			\
  const char* className##_Type = #className;
CLASSTYPE( CLASSNAME )

#define PARENT( classPath )							\
    INCLUDEPARENT( classPath )
#define __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
   ((className*)self)->_##methodName = _##className##_##methodName;
#define __VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
    __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#define __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
   ((className*)self)->_##methodName = NULL;
#define __VOIDABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
  __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )
#define OVERRIDE( methodName, type, argTypes, argNames )			\
  VIRTUALMETHOD( methodName, type, argTypes, argNames )
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )			\
  VOIDVIRTUALMETHOD( methodName, type, argTypes, argNames )
#define __METHOD( className, methodName, type, argTypes, argNames )
#define MEMBER( type, name )
#define _CLASSINIT( className )							\
   void className##_Init( void* self ) {					\
      assert( self );								\
      ((className*)self)->type = (char*)className##_Type;			\
      ((className*)self)->newFunc = (NewFunc*)className##_New;			\
      ((className*)self)->initFunc = className##_Init;
#define CLASSINIT( className )							\
   _CLASSINIT( className )
#define _CLASSCONSTRUCT( className )						\
  className##_Construct( self );
#define CLASSCONSTRUCT( className )						\
  _CLASSCONSTRUCT( className )
CLASSINIT( CLASSNAME )
#include INCLUDECLASS( CLASSNAME )
CLASSCONSTRUCT( CLASSNAME )
}

#define _CLASSNEW( className )							\
  className* className##_New() {						\
    className* self = (className*)MemAlloc( className, className##_Type );	\
    className##_Init( self );							\
    return self;								\
  }
#define CLASSNEW( className )							\
  _CLASSNEW( className )
CLASSNEW( CLASSNAME )

#undef PARENT
#define PARENT( classPath )							\
  INCLUDEFILE( StGermain/Base/Foundation, NoClass.h )
#undef __VIRTUALMETHOD
#define __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
  type className##_##methodName argTypes {					\
    assert( self );								\
    assert( ((className*)self)->_##methodName );				\
    return ((className*)self)->_##methodName argNames;				\
  }
#undef __VOIDVIRTUALMETHOD
#define __VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
  type className##_##methodName argTypes {					\
    assert( self );								\
    assert( ((className*)self)->_##methodName );				\
    ((className*)self)->_##methodName argNames;					\
  }
#undef __ABSTRACTMETHOD
#define __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
  __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#undef __VOIDABSTRACTMETHOD
#define __VOIDABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
  __VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#undef OVERRIDE
#define OVERRIDE( methodName, type, argTypes, argNames )			\
  VIRTUALMETHOD( methodName, type, argTypes, argNames )
#undef VOIDOVERRIDE
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )			\
   VOIDVIRTUALMETHOD( methodName, type, argTypes, argNames )
#define __METHOD( className, methodName, type, argTypes, argNames )
#undef MEMBER
#define MEMBER( type, name )
#include INCLUDECLASS( CLASSNAME )

#include "ClassPost.h"
