#define INCLUDEPATH( path )			\
  #path
#define INCLUDEFILE( dir, file )		\
  INCLUDEPATH( dir/file )
#define _INCLUDECLASS( dir, className )		\
  INCLUDEFILE( dir, className.def )
#define INCLUDECLASS( className )		\
  _INCLUDECLASS( CURRENTDIR, className )
#define INCLUDEHDR( classPath )			\
  INCLUDEPATH( classPath.h )
#define INCLUDEPARENT( classPath )		\
  INCLUDEPATH( classPath.def )

#define CLASSTYPE( className )	\
  _CLASSTYPE( className )

#define _VIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
  __VIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#define VIRTUALMETHOD( methodName, type, argTypes, argNames )			\
  _VIRTUALMETHOD( CLASSNAME, methodName, type, argTypes, argNames )
#define _VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )	\
  __VOIDVIRTUALMETHOD( className, methodName, type, argTypes, argNames )
#define VOIDVIRTUALMETHOD( methodName, type, argTypes, argNames )		\
  _VOIDVIRTUALMETHOD( CLASSNAME, methodName, type, argTypes, argNames )
#define _ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
  __ABSTRACTMETHOD( className, methodName, type, argTypes, argNames )
#define ABSTRACTMETHOD( methodName, type, argTypes, argNames )			\
  _ABSTRACTMETHOD( CLASSNAME, methodName, type, argTypes, argNames )
#define _VOIDABSTRACTMETHOD( className, methodName, type, argTypes, argNames )	\
  __VOIDABSTRACTMETHOD( className, methodName, type, argTypes, argNames )
#define VOIDABSTRACTMETHOD( methodName, type, argTypes, argNames )		\
  _VOIDABSTRACTMETHOD( CLASSNAME, methodName, type, argTypes, argNames )
#define _METHOD( className, methodName, type, argTypes, argNames )		\
  __METHOD( className, methodName, type, argTypes, argNames )
#define METHOD( methodName, type, argTypes, argNames )				\
  _METHOD( CLASSNAME, methodName, type, argTypes, argNames )
