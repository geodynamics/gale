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

#define CLASSTYPE( pre, className )		\
   _CLASSTYPE( pre, className )

#define _VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define VIRTUALMETHOD( methodName, type, argTypes, argNames )		\
   _VIRTUALMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define _VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VOIDVIRTUALMETHOD( pre, className, methodName, type, argTypes, argNames )
#define VOIDVIRTUALMETHOD( methodName, type, argTypes, argNames )	\
   _VOIDVIRTUALMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define _ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __ABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames )
#define ABSTRACTMETHOD( methodName, type, argTypes, argNames )		\
   _ABSTRACTMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define _VOIDABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames ) \
   __VOIDABSTRACTMETHOD( pre, className, methodName, type, argTypes, argNames )
#define VOIDABSTRACTMETHOD( methodName, type, argTypes, argNames )	\
   _VOIDABSTRACTMETHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
#define _METHOD( pre, className, methodName, type, argTypes, argNames )	\
   __METHOD( pre, className, methodName, type, argTypes, argNames )
#define METHOD( methodName, type, argTypes, argNames )			\
   _METHOD( PREFIX, CLASSNAME, methodName, type, argTypes, argNames )
