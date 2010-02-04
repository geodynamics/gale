Run all scripts from top level directory in stgUnderworld.

Run ./script/macroanalyze/createtables.sh first to create the necessary tables.
This creates 4 txt files.
defargs.txt  <- contains all macro DEFARG definitions.
parentchildhashtable.txt  <- contains all parent child _X_New function relations
proto.txt  <- contains all parent child _X_New function prototype arguments
structs.txt  <- contains all struct definitions

The scripts than can be directly called are:
-------------------------------------------------------------------------------------------------------------------
lineage.pl : uses parentchildhashtable.txt
Usage:
./script/macroanalyze/lineage _X_New
e.g.
./script/macroanalyze/lineage _PCDVC_New
will produce the hierarchy of functions
 _PCDVC_New => _DVCWeights_New => _WeightsCalculator_New => _Stg_Component_New => _Stg_Object_New => _Stg_Class_New

-------------------------------------------------------------------------------------------------------------------
readdefarg.pl : uses defargs.txt
Usage:
./script/macroanalyze/readdefarg.pl _X_DEFARGS
or
./script/macroanalyze/readdefarg.pl _X_New
e.g.
./script/macroanalyze/readdefarg.pl _PCDVC_New (or PCDVC_DEFARGS)
will give a list of the hierarchy of macros followed by the expanded form of the macro,
(with alternating colors of sets of args corresponding to additional args given by each parent macro)

PCDVC_DEFARGS
DVCWEIGHTS_DEFARGS
WEIGHTSCALCULATOR_DEFARGS
STG_COMPONENT_DEFARGS
STG_OBJECT_DEFARGS
STG_CLASS_DEFARGS
SizeT _sizeOfSelf,
Type type,
Stg_Class_DeleteFunction* _delete,
Stg_Class_PrintFunction* _print,
Stg_Class_CopyFunction* _copy

Name name,
AllocationType nameAllocationType 

Stg_Component_DefaultConstructorFunction* _defaultConstructor,
Stg_Component_ConstructFunction* _construct,
Stg_Component_BuildFunction* _build,
Stg_Component_InitialiseFunction* _initialise,
Stg_Component_ExecuteFunction* _execute,
Stg_Component_DestroyFunction* _destroy 

WeightsCalculator_CalculateFunction* _calculate 

-------------------------------------------------------------------------------------------------------------------
readstruct.pl : uses structs.txt
Usage:
./script/macroanalyze/readstruct.pl __X
e.g.
./script/macroanalyze/readstruct.pl __PCDVC
will  give a list of the hierarchy of structs followed by the expanded form of the struct,
(with alternating colors of sets of terms corresponding to additional terms given by each parent struct)

__PCDVC
__DVCWeights
__WeightsCalculator
__Stg_Component
__Stg_Object
__Stg_Class

SizeT _sizeOfSelf;
Type type;
Stg_Class_DeleteFunction* _delete;
Stg_Class_PrintFunction* _print;
Stg_Class_CopyFunction* _copy
Name name;
AllocationType nameAllocationType;
Stg_Component_DefaultConstructorFunction* _defaultConstructor;
Stg_Component_ConstructFunction* _construct;
Stg_Component_BuildFunction* _build;
Stg_Component_InitialiseFunction* _initialise;
Stg_Component_ExecuteFunction* _execute;
Stg_Component_DestroyFunction* _destroy;
Bool isConstructed;
Bool isBuilt;
Bool isInitialised;
Bool hasExecuted;
Bool isDestroyed;
Type constructType;
Type buildType;
Type initialiseType;
Type executeType;
Type destroyType;
FiniteElementContext* context;
WeightsCalculator_CalculateFunction* _calculate;
double cellLocalVolume;
Dimension_Index dim;
int resX;
int resY;
int resZ;
MaterialPointsSwarm* materialPointsSwarm;
double upperT;
double lowerT;
Bool splitInInterfaceCells;
Bool deleteInInterfaceCells;
int maxDeletions;
int maxSplits;
Bool Inflow;
double CentPosRatio;
int ParticlesPerCell;
double Threshold;
int maxDeletions_orig;
int maxSplits_orig;
Bool Inflow_orig;
Bool splitInInterfaceCells_orig;
Bool deleteInInterfaceCells_orig;

-------------------------------------------------------------------------------------------------------------------