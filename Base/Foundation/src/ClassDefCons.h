#define PARENT( classPath )							\
    INCLUDEPARENT( classPath )
#define __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   ((pre##className*)self)->_##methodName = _##pre##className##_##methodName;
#define __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames )	\
   ((pre##className*)self)->_##methodName = NULL;
#define __VOIDABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames )
#define OVERRIDE( methodName, type, argTypes, argNames )	\
   _VIRTUALMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define VOIDOVERRIDE( methodName, type, argTypes, argNames )	\
   _VOIDVIRTUALMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define __METHOD( pre, className, methodName, type, argTypes, argNames )
#define MEMBER( type, name )
#define _CLASSCONSTRUCT( pre, className )				\
   void pre##className##_Construct( void* self ) {				\
      assert( self );								\
      ((pre##className*)self)->type = (char*)pre##className##_Type;		\
      ((pre##className*)self)->newFunc = (NewFunc*)pre##className##_New;	\
      ((pre##className*)self)->constructFunc = pre##className##_Construct;
#define CLASSCONSTRUCT( pre, className )		\
   _CLASSCONSTRUCT( pre, className )
#define _CLASSINIT( pre, className )		\
  pre##className##_Init( self );
#define CLASSINIT( pre, className )		\
   _CLASSINIT( pre, className )

CLASSCONSTRUCT( PREFIX, CLASSNAME )
#include INCLUDECLASS( CLASSNAME )
CLASSINIT( PREFIX, CLASSNAME )
}

#undef PARENT
#undef __VIRTUALMETHOD
#undef __VOIDVIRTUALMETHOD
#undef __ABSTRACTMETHOD
#undef __VOIDABSTRACTMETHOD
#undef OVERRIDE
#undef VOIDOVERRIDE
#undef __METHOD
#undef MEMBER
#undef _CLASSCONSTRUCT
#undef CLASSCONSTRUCT
#undef _CLASSINIT
#undef CLASSINIT
