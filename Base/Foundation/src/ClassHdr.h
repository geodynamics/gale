#ifndef CLASSNAME
#error No class name given.
#endif

#include "ClassPre.h"

/*
#define PARENT( classPath )							\
    INCLUDEHDR( classPath )
#define __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#define __VOIDVIRTUALMETHOD( className, methodName, argTypes, argNames )
#define __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )
#define __VOIDABSTRACTMETHOD( className, methodName, argTypes, argNames )
#define __OVERRIDE( className, methodName, type, argTypes, argNames )
#define __METHOD( className, methodName, type, argTypes, argNames )
#define MEMBER( type, name )
#include INCLUDECLASS( CLASSNAME )
*/

#define _CLASSTYPE( className )		\
  extern const char* className##_Type;
CLASSTYPE( CLASSNAME )

#define PARENT( classPath )							\
  INCLUDEFILE( StGermain/Base/Foundation, NoClass.h )
#define __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
  typedef type (className##_##methodName##Func)argTypes;
#define __VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
    __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#define __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
    __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#define __VOIDABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
    __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#define OVERRIDE( methodName, type, argTypes, argNames )
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )
#define __METHOD( className, methodName, type, argTypes, argNames )
#define MEMBER( type, name )
#include INCLUDECLASS( CLASSNAME )

#undef PARENT
#define PARENT( classPath )							\
  INCLUDEPARENT( classPath )
#undef __VIRTUALMETHOD
#define __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
  className##_##methodName##Func*	_##methodName;
#undef __VOIDVIRTUALMETHOD
#define __VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
    __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#undef __ABSTRACTMETHOD
#define __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
    __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#undef __VOIDABSTRACTMETHOD
#define __VOIDABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
    __VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#undef MEMBER
#define MEMBER( type, name )							\
  type	name;
struct CLASSNAME {
#include INCLUDECLASS( CLASSNAME )
};

#define _CLASSNEW( className )			\
  className* className##_New();
#define CLASSNEW( className )			\
  _CLASSNEW( className )
#define _CLASSINIT( className )			\
  void className##_Init( void* self );
#define CLASSINIT( className )			\
  _CLASSINIT( className )
CLASSNEW( CLASSNAME )
CLASSINIT( CLASSNAME )

#undef PARENT
#define PARENT( classPath )							\
  INCLUDEFILE( StGermain/Base/Foundation, NoClass.h )
#undef __VIRTUALMETHOD
#define __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
    type _##className##_##methodName argTypes;					\
    type className##_##methodName argTypes;
#undef __VOIDVIRTUALMETHOD
#define __VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
    __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#undef __ABSTRACTMETHOD
#define __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
  type className##_##methodName argTypes;
#undef __VOIDABSTRACTMETHOD
#define __VOIDABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
    __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )
#undef OVERRIDE
#define OVERRIDE( methodName, type, argTypes, argNames )			\
  VIRTUALMETHOD( methodName, type, argTypes, argNames )
#undef VOIDOVERRIDE
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )			\
  VOIDVIRTUALMETHOD( methodName, type, argTypes, argNames )
#undef __METHOD
#define __METHOD( className, methodName, type, argTypes, argNames )		\
  type className##_##methodName argTypes;
#undef MEMBER
#define MEMBER( type, name )
#include INCLUDECLASS( CLASSNAME )

#include "ClassPost.h"
