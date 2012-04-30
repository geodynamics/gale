/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** $Id: BuoyancyForceTerm.c 518 2007-10-11 08:07:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>

#include "types.h"
#include "BuoyancyForceTerm.h"
#include "HydrostaticTerm.h"
#include "MaterialSwarmVariable.h"
#include <list>
#include <vector>
#include <string>

#include <assert.h>
#include <string.h>
#include <stddef.h>
#include <sstream>

static std::vector<std::string> equations[2];

/* Textual name of this class */
const Type BuoyancyForceTerm_Type = "BuoyancyForceTerm";

BuoyancyForceTerm* BuoyancyForceTerm_New( 
	Name							name,
	FiniteElementContext*	context,
	ForceVector*				forceVector,
	Swarm*						integrationSwarm,
	FeVariable*					temperatureField,
	FeVariable*					pressureField,
	double						gravity,
        double*                                         gHat,
	Bool							adjust,
	Bool							damping,
	Materials_Register*		materials_Register,
        HydrostaticTerm*                                hydrostaticTerm )
{
	BuoyancyForceTerm* self = (BuoyancyForceTerm*) _BuoyancyForceTerm_DefaultNew( name );

	self->isConstructed = True;
	_ForceTerm_Init( self, context, forceVector, integrationSwarm, NULL );
	_BuoyancyForceTerm_Init(self,temperatureField,pressureField,gravity,
                                gHat,adjust,damping,materials_Register,
                                hydrostaticTerm);

	return self;
}

/* Creation implementation / Virtual constructor */
BuoyancyForceTerm* _BuoyancyForceTerm_New(  BUOYANCYFORCETERM_DEFARGS  )
{
	BuoyancyForceTerm* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(BuoyancyForceTerm) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (BuoyancyForceTerm*) _ForceTerm_New(  FORCETERM_PASSARGS  );
	
	/* Virtual info */
	self->_calcGravity = _calcGravity;
	
	return self;
}

void _BuoyancyForceTerm_Init( 
	void*						forceTerm, 
	FeVariable*				temperatureField,
	FeVariable*				pressureField,
	double					gravity,
        double*                                 gHat,
	Bool						adjust,
	Bool						damping,
	Materials_Register*	materials_Register,
        HydrostaticTerm*                        hydrostaticTerm )
{
	BuoyancyForceTerm* self = (BuoyancyForceTerm*)forceTerm;

	self->temperatureField    = temperatureField;
	self->pressureField       = pressureField;
	self->gravity             = gravity;
	self->gHat		  = gHat;
	self->adjust              = adjust;
	self->damping             = damping;
	self->materials_Register  = materials_Register;
        self->hydrostaticTerm     = hydrostaticTerm;
}

void _BuoyancyForceTerm_Delete( void* forceTerm ) {
	BuoyancyForceTerm* self = (BuoyancyForceTerm*)forceTerm;

	_ForceTerm_Delete( self );
}

void _BuoyancyForceTerm_Print( void* forceTerm, Stream* stream ) {
	BuoyancyForceTerm* self = (BuoyancyForceTerm*)forceTerm;
	
	_ForceTerm_Print( self, stream );

	/* General info */
	Journal_PrintPointer( stream, self->temperatureField );
	Journal_PrintPointer( stream, self->pressureField );
	Journal_PrintDouble( stream, self->gravity );
}

void* _BuoyancyForceTerm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(BuoyancyForceTerm);
	Type                                                      type = BuoyancyForceTerm_Type;
	Stg_Class_DeleteFunction*                              _delete = _BuoyancyForceTerm_Delete;
	Stg_Class_PrintFunction*                                _print = _BuoyancyForceTerm_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _BuoyancyForceTerm_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _BuoyancyForceTerm_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _BuoyancyForceTerm_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _BuoyancyForceTerm_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _BuoyancyForceTerm_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _BuoyancyForceTerm_Destroy;
	ForceTerm_AssembleElementFunction*            _assembleElement = _BuoyancyForceTerm_AssembleElement;
	BuoyancyForceTerm_CalcGravityFunction*            _calcGravity = _BuoyancyForceTerm_CalcGravity;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*)_BuoyancyForceTerm_New(  BUOYANCYFORCETERM_PASSARGS  );
}

void _BuoyancyForceTerm_AssignFromXML(void* forceTerm,Stg_ComponentFactory* cf,
                                      void* data)
{
  BuoyancyForceTerm*		self = (BuoyancyForceTerm*)forceTerm;
  Dictionary*					dict;
  Dictionary_Entry_Value*	tmp;
  char*							rootKey;
  FeVariable*					temperatureField;
  FeVariable*					pressureField;
  double						gravity;
  Bool							adjust;
  Bool							damping;
  Materials_Register*		materials_Register;
  unsigned						nDims;
  Dictionary_Entry_Value*	direcList;
  double*						gHat;
  unsigned						d_i;
  PICelleratorContext*		context;
  HydrostaticTerm*                                    hydrostaticTerm;

  /* Construct Parent */
  _ForceTerm_AssignFromXML( self, cf, data );

  dict=Dictionary_Entry_Value_AsDictionary(Dictionary_Get
                                           (cf->componentDict,
                                            (Dictionary_Entry_Key)self->name));
                                                          
  temperatureField=
    Stg_ComponentFactory_ConstructByKey(cf,self->name,
                                        (Dictionary_Entry_Key)"TemperatureField",
                                        FeVariable, False, data  ) ;
  pressureField=
    Stg_ComponentFactory_ConstructByKey(cf,self->name,
                                        (Dictionary_Entry_Key)"PressureField",
                                        FeVariable,False,data);
  gravity=
    Stg_ComponentFactory_GetDouble(cf,self->name,
                                   (Dictionary_Entry_Key)"gravity",0.0);
  adjust=Stg_ComponentFactory_GetBool(cf,self->name,
                                      (Dictionary_Entry_Key)"adjust",False);
  damping=Stg_ComponentFactory_GetBool(cf,self->name,
                                       (Dictionary_Entry_Key)"damping",True);

  direcList=Dictionary_Get(dict,(Dictionary_Entry_Key)"gravityDirection");

  if(direcList)
    {
      nDims = Dictionary_Entry_Value_GetCount( direcList  );
      gHat = AllocArray( double, nDims );

      for( d_i = 0; d_i < nDims; d_i++ )
        {
          tmp = Dictionary_Entry_Value_GetElement( direcList, d_i );
          rootKey = Dictionary_Entry_Value_AsString( tmp );

          if( !Stg_StringIsNumeric( (char *)rootKey )  )
            tmp = Dictionary_Get( cf->rootDict, (Dictionary_Entry_Key)rootKey );
          gHat[d_i] = Dictionary_Entry_Value_AsDouble( tmp );
        }
      if( nDims == 2  )
        Vec_Norm2D( gHat, gHat );
      else
        Vec_Norm3D( gHat, gHat );
    }
  else
    {
      gHat = AllocArray( double, 3 );
      gHat[0]=0;
      gHat[1]=-1;
      gHat[2]=0;
    }

  context = (PICelleratorContext*)self->context;
  assert( Stg_CheckType( context, PICelleratorContext ) );
  materials_Register = context->materials_Register;
  assert( materials_Register );

  hydrostaticTerm = Stg_ComponentFactory_ConstructByKey(cf, self->name,
                                                        "HydrostaticTerm",
                                                        HydrostaticTerm, False,
                                                        data);

  _BuoyancyForceTerm_Init(self,temperatureField,pressureField,gravity,gHat,
                          adjust,damping,materials_Register,hydrostaticTerm);
}

void _BuoyancyForceTerm_Build(void* forceTerm, void* data)
{
  BuoyancyForceTerm* self = (BuoyancyForceTerm*)forceTerm;
  BuoyancyForceTerm_MaterialExt*   materialExt;
  Material_Index                   material_I;
  Material*                        material;
  Materials_Register* materials_Register = self->materials_Register;
  IntegrationPointsSwarm* swarm=(IntegrationPointsSwarm*)self->integrationSwarm;
  MaterialPointsSwarm**            materialSwarms;
  Index                            materialSwarm_I;
  char*                            name;
  Stg_ComponentFactory*            cf;

  cf = self->context->CF;

  _ForceTerm_Build( self, data );

  if ( self->temperatureField )
    Stg_Component_Build( self->temperatureField, data, False );
  if ( self->pressureField )
    Stg_Component_Build( self->pressureField, data, False );

  /* Sort out material extension stuff */
  self->materialExtHandle=
    Materials_Register_AddMaterialExtension(self->materials_Register, 
                                            self->type,
                                            sizeof(BuoyancyForceTerm_MaterialExt));

  Journal_Firewall(equations[0].empty() && equations[1].empty(),
                   Journal_MyStream(Error_Type,self),
                   "You may only create one instance of BuoyancyForceTerm");
          
  for(material_I=0;
      material_I<Materials_Register_GetCount(materials_Register);
      material_I++)
    {
      material=Materials_Register_GetByIndex(materials_Register,material_I);
      materialExt=(BuoyancyForceTerm_MaterialExt*)
        ExtensionManager_GetFunc(material->extensionMgr,material,
                                 self->materialExtHandle);

      /* Set the density */
      equations[0].push_back(Stg_ComponentFactory_GetString
                             (cf,material->name,
                              (Dictionary_Entry_Key)"densityEquation",
                              ""));
                   
      materialExt->density =
        Stg_ComponentFactory_GetDouble
        (cf,material->name,(Dictionary_Entry_Key)"density",0.0);

      /* Set alpha */
      equations[1].push_back(Stg_ComponentFactory_GetString
                             (cf,material->name,
                              (Dictionary_Entry_Key)"alphaEquation",
                              ""));
                   
      materialExt->alpha =
        Stg_ComponentFactory_GetDouble
        (cf,material->name,(Dictionary_Entry_Key)"alpha",0.0);
    }
	
  /* Create Swarm Variables of each material swarm this ip swarm
     is mapped against */
  materialSwarms=IntegrationPointMapper_GetMaterialPointsSwarms
    (swarm->mapper, &(self->materialSwarmCount) );
  self->densitySwarmVariables=
    Memory_Alloc_Array(MaterialSwarmVariable*, self->materialSwarmCount,
                       "DensityVariables" );
  self->alphaSwarmVariables=
    Memory_Alloc_Array(MaterialSwarmVariable*, self->materialSwarmCount,
                       "AlphaVariables");
	
  for(materialSwarm_I=0; materialSwarm_I<self->materialSwarmCount;
      ++materialSwarm_I )
    {
      name = Stg_Object_AppendSuffix( materialSwarms[materialSwarm_I], (Name)"Density"  );
      self->densitySwarmVariables[materialSwarm_I]=
        MaterialSwarmVariable_New(name,
                                  (AbstractContext*)self->context,
                                  materialSwarms[materialSwarm_I], 
                                  1, 
                                  self->materials_Register, 
                                  self->materialExtHandle, 
                                  GetOffsetOfMember( *materialExt, density ) );
      Memory_Free( name );

      name=Stg_Object_AppendSuffix(materialSwarms[materialSwarm_I],
                                   (Name)"Alpha");
      self->alphaSwarmVariables[materialSwarm_I]=
        MaterialSwarmVariable_New(name,
                                  (AbstractContext*)self->context,
                                  materialSwarms[materialSwarm_I], 
                                  1, 
                                  self->materials_Register, 
                                  self->materialExtHandle, 
                                  GetOffsetOfMember( *materialExt, alpha ) );
      Memory_Free( name );

      /* Build new Swarm Variables */
      Stg_Component_Build(self->densitySwarmVariables[materialSwarm_I],
                          data, False);
      Stg_Component_Build(self->alphaSwarmVariables[materialSwarm_I],
                          data, False);
    }
}

void _BuoyancyForceTerm_Initialise( void* forceTerm, void* data )
{
  BuoyancyForceTerm* self=(BuoyancyForceTerm*)forceTerm;
  Index i;

  _ForceTerm_Initialise( self, data );

  if ( self->temperatureField )
    Stg_Component_Initialise( self->temperatureField, data, False );
  if ( self->pressureField )
    Stg_Component_Initialise( self->pressureField, data, False );
	
  for ( i = 0; i < self->materialSwarmCount; ++i ) {
    Stg_Component_Initialise( self->densitySwarmVariables[i], data, False );
    Stg_Component_Initialise( self->alphaSwarmVariables[i],   data, False );
  }
}

void _BuoyancyForceTerm_Execute( void* forceTerm, void* data )
{
  BuoyancyForceTerm* self = (BuoyancyForceTerm*)forceTerm;
  _ForceTerm_Execute( self, data );
}

void _BuoyancyForceTerm_Destroy( void* forceTerm, void* data )
{
  BuoyancyForceTerm* self = (BuoyancyForceTerm*)forceTerm;
  Index i;

  for ( i = 0; i < self->materialSwarmCount; ++i ) {
    _Stg_Component_Delete( self->densitySwarmVariables[i] );
    _Stg_Component_Delete( self->alphaSwarmVariables[i] );
  }

  FreeArray( self->gHat );

  Memory_Free( self->densitySwarmVariables );
  Memory_Free( self->alphaSwarmVariables );

  _ForceTerm_Destroy( forceTerm, data );
}

void _BuoyancyForceTerm_AssembleElement(void* forceTerm,
                                        ForceVector* forceVector,
                                        Element_LocalIndex lElement_I,
                                        double* elForceVec)
{
  BuoyancyForceTerm* self = (BuoyancyForceTerm*) forceTerm;
  Dimension_Index dim = forceVector->dim;
  IntegrationPointsSwarm* swarm=(IntegrationPointsSwarm*)self->integrationSwarm;
  FeMesh* mesh = forceVector->feVariable->feMesh;

  ElementType* elementType = FeMesh_GetElementType( mesh, lElement_I );
  Element_NodeIndex elementNodeCount = elementType->nodeCount;
  Dof_Index nodeDofCount = dim;
  Cell_Index cell_I=
    CellLayout_MapElementIdToCellId(swarm->cellLayout,lElement_I);
  Particle_InCellIndex cellParticleCount = swarm->cellParticleCountTbl[cell_I];

  /* adjust & adjustFactor -- 20060411 Alan
   *
   * The adjust decides whether an adjustFactor should be applied to
   * the resulting factor.  If on, the total weight of the particles
   * in the cell are scaled to the cell local volume.
   *
   * This is designed to be used when integrating with swarms which do
   * not cover the whole domain (ie - I used it to do dave.m's test of
   * 1 swarm for blob, 1 swarm for background)
   */ 
  double adjustFactor;
  if(self->adjust)
    {
      double totalWeight = 0.0;
      for(Particle_InCellIndex cParticle_I=0; cParticle_I<cellParticleCount;
          cParticle_I++)
        {
          IntegrationPoint* particle=
            (IntegrationPoint*)Swarm_ParticleInCellAt(swarm,cell_I,cParticle_I);
          totalWeight += particle->weight;
        }
      adjustFactor = swarm->weights->cellLocalVolume / totalWeight;
    }
  else
    {
      adjustFactor = 1.0;
    }

  /* Use an average density and alpha for OneToMany mapper */
  double density(0), alpha(0);
  bool one_to_many=Stg_Class_IsInstance(swarm->mapper,OneToManyMapper_Type);
  if(one_to_many)
    {
      BuoyancyForceTerm_average_density_alpha(self,swarm,mesh,lElement_I,
                                              &density,&alpha);
    }

  for(Particle_InCellIndex cParticle_I=0; cParticle_I<cellParticleCount;
      cParticle_I++)
    {
      IntegrationPoint* particle=
        (IntegrationPoint*)Swarm_ParticleInCellAt(swarm,cell_I,cParticle_I);
      double *xi=particle->xi;

      double detJac=ElementType_JacobianDeterminant(elementType,mesh,
                                                    lElement_I,xi,dim);
      double Ni[27];
      ElementType_EvaluateShapeFunctionsAt(elementType,xi,Ni);

      /* Get parameters */
      double temperature(0), pressure(0);
      if(self->temperatureField) 
        FeVariable_InterpolateFromMeshLocalCoord
          (self->temperatureField, mesh, lElement_I, xi, &temperature);

      if(self->pressureField)
        FeVariable_InterpolateFromMeshLocalCoord
          (self->pressureField, mesh, lElement_I, xi, &pressure);

      double background_density = 0.0;
      if(self->hydrostaticTerm)
        {
          Coord coord;
          FeMesh_CoordLocalToGlobal(mesh, cell_I, xi, coord);
          background_density=
            HydrostaticTerm_Density(self->hydrostaticTerm,coord);
          if(self->pressureField)
            pressure+=HydrostaticTerm_Pressure(self->hydrostaticTerm,
                                               coord);
        }

      if(!one_to_many)
        {
          /* Handle case where we are using gauss swarms with
             NearestNeighborMapper instead of a material
             swarm */
          IntegrationPointsSwarm* NNswarm(swarm);
          IntegrationPoint* NNparticle(particle);
          NearestNeighbor_Replace(&NNswarm,&NNparticle,lElement_I,dim);

          /* Get density and alpha */
          density=BuoyancyForceTerm_density(self,mesh,lElement_I,cell_I,
                                            NNswarm->mapper,NNparticle);
          alpha=BuoyancyForceTerm_alpha(self,mesh,lElement_I,cell_I,
                                        NNswarm->mapper,NNparticle);
        }
      /* Calculate Force */
      double gravity[dim];
      BuoyancyForceTerm_CalcGravity(self,mesh,lElement_I,xi,gravity);

      double rho=(density*(1.0-alpha*temperature) - background_density);
      double factor = detJac * particle->weight * adjustFactor * rho;

      /* Apply force in the correct direction */
      for(Node_ElementLocalIndex eNode_I=0; eNode_I<elementNodeCount; eNode_I++)
        {
          for(unsigned d_i=0; d_i<dim; d_i++)
            elForceVec[ eNode_I * nodeDofCount + d_i ] +=
              gravity[d_i] * factor * Ni[ eNode_I ] ;
        }
    }
}

/* The default implementation is for the gravity to be constant. */
void _BuoyancyForceTerm_CalcGravity(BuoyancyForceTerm *self, FeMesh *mesh,
                                    Element_DomainIndex dElement_I,
                                    double* xi, double *gravity)
{
  for(unsigned int d=0;d<mesh->generator->nDims;++d)
    gravity[d]=self->gHat[d];
}

void BuoyancyForceTerm_average_density_alpha(BuoyancyForceTerm *self,
                                             IntegrationPointsSwarm* input_swarm,
                                             FeMesh* mesh,
                                             const Element_LocalIndex &lElement_I,
                                             double *density, double *alpha)
{
  *density=*alpha=0;
  IntegrationPointsSwarm* swarm(input_swarm);
  bool one_to_many=Stg_Class_IsInstance(swarm->mapper,OneToManyMapper_Type);
  if(one_to_many)
    {
      swarm=((OneToManyMapper*)(swarm->mapper))->swarm;
    }
  int cell_I=CellLayout_MapElementIdToCellId(swarm->cellLayout,
                                             lElement_I);
  int num_particles=swarm->cellParticleCountTbl[cell_I];

  double weight(0);
  for(int ii=0;ii<num_particles;++ii)
    {
      IntegrationPoint *particle=
        (IntegrationPoint*)Swarm_ParticleInCellAt(swarm,cell_I,ii);
                                                  
      weight+=particle->weight;
      *density+=
        BuoyancyForceTerm_density(self,mesh,lElement_I,cell_I,
                                  swarm->mapper,particle)*particle->weight;

      *alpha+=BuoyancyForceTerm_alpha(self,mesh,lElement_I,cell_I,
                                      swarm->mapper,particle)*particle->weight;
        
    }
  *density/=weight;
  *alpha/=weight;
}

double BuoyancyForceTerm_density_alpha(BuoyancyForceTerm *self,
                                       FeMesh* mesh,
                                       const Element_LocalIndex &lElement_I,
                                       const Cell_Index &cell_I,
                                       IntegrationPointMapper* mapper,
                                       IntegrationPoint* particle,
                                       const int &type)
{
  unsigned int material_index=
    IntegrationPointMapper_GetMaterialIndexOn(mapper,particle);

  Journal_Firewall(material_index<equations[type].size(),
                   Journal_MyStream(Error_Type,self),
                   "Bad material index %d",material_index);

  std::string equation = equations[type][material_index];
  double result;
  if(!equation.empty())
    {
      std::stringstream ss;

      double temperature(0), pressure(0);
      if(self->temperatureField) 
        FeVariable_InterpolateFromMeshLocalCoord
          (self->temperatureField, mesh, lElement_I, particle->xi, &temperature);

      if(self->pressureField)
        FeVariable_InterpolateFromMeshLocalCoord
          (self->pressureField, mesh, lElement_I, particle->xi, &pressure);

      Coord coord;
      FeMesh_CoordLocalToGlobal(mesh, cell_I, particle->xi, coord);

      ss << "T=" << temperature
         << "; p=" << pressure
         << "; " << equation;
      result=Equation_eval(coord,
                           (DomainContext *)(self->context),
                           ss.str());
    }
  else
    {
      if(type==0)
        {
          result=IntegrationPointMapper_GetDoubleFromMaterial
            (mapper, particle, self->materialExtHandle,
             offsetof(BuoyancyForceTerm_MaterialExt, density));
        }
      else
        {
          result=IntegrationPointMapper_GetDoubleFromMaterial
            (mapper, particle, self->materialExtHandle,
             offsetof(BuoyancyForceTerm_MaterialExt, alpha));
        }
    }
  return result;
}
